"""Score a SanntiS bundle on the nucleotide test sets (04_testsets.sh) with the SanntiS metric code: Test A / B on
the isolate and shuffled contigs at the bundle's t* (metrics with 95 % CIs over source contigs, strata by class /
lineage / novelty, radar, PR curves and AUPRC, NPC classification on true loci and predicted regions), and the
evaluation tables objective 03 reads (``test_bgcs.parquet`` with tokens and novelty, ``labels/``).

    python score_testsets.py --testsets DIR --test-bgcs PARQUET --tokens PARQUET --model BUNDLE
        --train-entities PARQUET --labels-loci PARQUET --out DIR [--expected DIR] [--n-boot 1000]

``--train-entities``: the training run's ``entities/entities.parquet`` (its positives' Pfam sets define novelty:
max Dice < 0.5 novel, < 0.9 related, else near-duplicate). ``--labels-loci``: a ``collect/loci.parquet`` holding
the test BGCs' NPClassifier labels (objective 01's run: MIBiG 5.0rc1 active + pending). ``--expected``: the
release-run tables, reported side by side as a reference (other backgrounds, not a target).
Writes <out>/{test_bgcs.parquet, labels/, score/{A,B}-{isolate,shuffled}.pkl, score/test_*.tsv, score/summary.tsv,
score/reference.tsv, summary.json}.
"""

import argparse
import json
import logging
import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sanntis.config import load_recipe
from sanntis.model.bundle import SanntisModel
from sanntis.train import metrics as M
from sanntis.train.labels import superclass_parents

from sv1 import expected as E
from sv1.archives import sha256
from sv1.novelty import max_dice, novelty_bin
from sv1.testscore import evaluate_dataset, npc_evaluate, predict_contigs
from sv1.testsets import contigs, load_cds, sets_of

log = logging.getLogger("score_testsets")
HEADLINE = ("precision_t", "recall_macro_t", "f1_macro_t", "bgc_precision", "bgc_macro_recall")


def lists(x) -> list:
    return [] if x is None else [lists(v) if isinstance(v, np.ndarray | list) else v for v in x]


def labels_of(loci: Path, out: Path, status: pd.Series) -> tuple[pd.DataFrame, dict]:
    """NPC labels of the test BGCs from a ``collect/loci.parquet`` (as ``sanntis train`` labels them), written in
    the ``labels/`` layout of the ``v3active_eval`` deposit; ``status`` = MIBiG entry status per test entity_id."""
    ids = list(status.index)
    df = pd.read_parquet(loci)
    df = df[df["entity_id"].isin(ids)]
    missing = sorted(set(ids) - set(df["entity_id"]))
    if missing:
        raise SystemExit(f"{len(missing)} test BGCs without a row in {loci}: {missing[:5]}")
    lab = pd.DataFrame({"accession": df["entity_id"].to_numpy(),
                        "curation_status": df["entity_id"].map(status).to_numpy(),
                        "npc_status": df["npc_status"].fillna("unlabelled").to_numpy(),
                        "pathway": [lists(x) for x in df["npc_pathway"]],
                        "superclass": [lists(x) for x in df["npc_superclass"]]}).sort_values("accession")
    parents = {k: sorted(v) for k, v in superclass_parents().items()}
    d = out / "labels"
    d.mkdir(parents=True, exist_ok=True)
    lab.to_parquet(d / "npc_labels.parquet", index=False)
    (d / "npc_superclass_parents.json").write_text(json.dumps(parents, indent=1) + "\n")
    (d / "manifest.json").write_text(json.dumps({
        "source": {"path": str(loci), "sha256": sha256(loci)}, "n": len(lab),
        "npc_status": lab["npc_status"].value_counts().to_dict()}, indent=1) + "\n")
    return lab.set_index("accession"), parents


def test_table(tb: pd.DataFrame, testsets: Path, tokens: dict, pfam: dict, train_entities: Path) -> pd.DataFrame:
    """Tokens per CDS (gene order of the archive locus, from the isolate contigs), Pfam set, novelty vs training."""
    toks, pf = {}, {}
    for _, c in load_cds(testsets, "isolate").values():
        for e, g in c[c["source_kind"] == "positive"].groupby("entity_id"):
            g = g.sort_values("start")
            toks[e] = [tokens.get(m, []) for m in g["md5"]]
            pf[e] = sorted(set().union(*(set(pfam.get(m, [])) for m in g["md5"])))
    tb = tb.copy()
    tb["cds_tokens"] = [toks[e] for e in tb["entity_id"]]
    tb["pfam"] = [pf[e] for e in tb["entity_id"]]
    tr = pd.read_parquet(train_entities, columns=["entity_id", "kind", "pfam"])
    tr = tr[tr["kind"] == "positive"]
    d = max_dice([set(x) for x in tb["pfam"]], [set(lists(x)) for x in tr["pfam"]])
    tb["max_dice"] = d
    tb["novelty"] = [novelty_bin(x) for x in d]
    tb["classes"] = [lists(x) for x in tb["classes"]]
    no_tok = int(sum(not any(t) for t in tb["cds_tokens"]))
    log.info("test BGCs: %d, novelty %s, without any token %d", len(tb), tb["novelty"].value_counts().to_dict(),
             no_tok)
    return tb


def reference(score: Path, exp: Path | None) -> pd.DataFrame | None:
    if not exp or not (exp / "test_metrics.tsv").exists():
        return None
    new = pd.read_csv(score / "test_metrics.tsv", sep="\t")
    old = pd.read_csv(exp / "test_metrics.tsv", sep="\t")
    key = ["set", "scope", "metric"]
    m = new.merge(old, on=key, suffixes=("", "_release"), how="left")
    m = m[(m["scope"] == "all") & m["metric"].isin(HEADLINE)]
    return m[[*key, "value", "ci_lo", "ci_hi", "value_release", "ci_lo_release", "ci_hi_release"]]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    p.add_argument("--testsets", type=Path, required=True)
    p.add_argument("--test-bgcs", type=Path, required=True)
    p.add_argument("--tokens", type=Path, required=True)
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--train-entities", type=Path, required=True)
    p.add_argument("--labels-loci", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--expected", type=Path, default=None)
    p.add_argument("--n-boot", type=int, default=1000)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    out, score = a.out, a.out / "score"
    score.mkdir(parents=True, exist_ok=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = SanntisModel.from_pretrained(a.model)
    t_star = float(model.config.decode.threshold)
    grid = sorted(set(M.threshold_grid(load_recipe().operating_point.grid)) | {t_star})
    tok = pd.read_parquet(a.tokens)
    tokens = dict(zip(tok["md5"], map(lists, tok["tokens"]), strict=True))
    pfam = dict(zip(tok["md5"], map(lists, tok["pfam"]), strict=True))

    tb = test_table(pd.read_parquet(a.test_bgcs), a.testsets, tokens, pfam, a.train_entities)
    tb.to_parquet(out / "test_bgcs.parquet", index=False)
    meta = tb.set_index("entity_id").assign(classes=lambda d: d["classes"].map(tuple))
    pools = {"A": set(tb.loc[tb["test_a"], "entity_id"]), "B": set(tb.loc[tb["test_b"], "entity_id"])}
    labels, parents = labels_of(a.labels_loci, out, meta.loc[sorted(pools["B"]), "status"])

    summary = {"model": str(a.model), "t_star": t_star, "test_A": len(pools["A"]), "test_B": len(pools["B"])}
    rows, extra = [], {}
    for kind in sets_of(a.testsets):
        cds = load_cds(a.testsets, kind)
        cs = contigs(cds, tokens)
        preds = predict_contigs(model, cs, grid)
        for test, ids in pools.items():
            name = f"{test}-{kind}"
            res = evaluate_dataset(preds, t_star, grid, meta, ids, n_boot=a.n_boot, shuffled=kind == "shuffled")
            if model.npc is not None:
                res["npc"], res["npc_per_pathway"] = npc_evaluate(model, preds, t_star, tb, labels, parents, ids)
            res["n_contigs"] = len(cs)
            (score / f"{name}.pkl").write_bytes(pickle.dumps(res))
            summary[name] = {k: round(float(res["metrics"][k]), 4) for k in HEADLINE}
            log.info("%s: %s", name, summary[name])
            rows += E.metric_rows(res, set=name)
            rows.append({"set": name, "scope": "all", "metric": "n_contigs", "value": res["n_contigs"]})
            for k, df in E.frames(res, set=name).items():
                extra.setdefault(k, []).append(df)
    pd.DataFrame(rows).to_csv(score / "test_metrics.tsv", sep="\t", index=False, float_format="%.6g")
    for k, dfs in extra.items():
        pd.concat(dfs).to_csv(score / f"test_{k}.tsv", sep="\t", index=False, float_format="%.6g")
    pd.DataFrame([{"set": s, **v} for s, v in summary.items() if isinstance(v, dict)]).to_csv(
        score / "summary.tsv", sep="\t", index=False)
    ref = reference(score, a.expected)
    if ref is not None:
        ref.to_csv(score / "reference.tsv", sep="\t", index=False, float_format="%.4g")
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    log.info("wrote %s", score)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
