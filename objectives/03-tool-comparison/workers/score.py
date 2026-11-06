"""Score SanntiS 1.0 and the other tools on the nucleotide test sets, every method on the same contigs and CDS,
with the SanntiS metric code (`sanntis.train.metrics`).

    python score.py --testsets DIR --tools DIR --model BUNDLE --test-bgcs PARQUET --tokens PARQUET --out DIR

``--tools`` (repeatable, optional) holds one directory per tool run (``<run>/<set>_<track>/``, the layout of
``02_tools.sh``); a tool whose directory is missing is left out. antiSMASH is lineage-matched as a user would run
it: ``--taxon bacteria`` on prokaryotic contigs, ``--taxon fungi`` (fungiSMASH) on eukaryotic ones, plantiSMASH 1.0
for regions over a plant BGC. Writes ``results.pkl`` (per test set and tool: metrics, CIs, strata), ``summary.tsv``
(one row per set and tool) and the flat tables of ``sv1.expected`` (metrics with CIs, strata, radar, PR curves,
AUPRC).
"""

import argparse
import logging
import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from parsers import prophet_regions, sanntis09_regions
from sanntis.config import load_recipe
from sanntis.model.bundle import SanntisModel
from sanntis.train import metrics as M
from testsets import KIND, SETS, TRACKS, contigs, load_cds

from sv1 import expected
from sv1 import metrics as SM
from sv1 import tools as T
from sv1.testscore import evaluate_dataset, predict_contigs, scored

LEVELS = ("relaxed", "strict", "loose")  # report order: relaxed is the reference row
log = logging.getLogger("score")


# runs scored as-is (no lineage matching), in report order: (run directory, report label, files, parser)
PLAIN = [("gecco_0.9.5", "GECCO 0.9.5", "*.clusters.tsv", T.gecco_regions),
         ("gecco_0.11.0", "GECCO 0.11.0", "*.clusters.tsv", T.gecco_regions),
         ("sanntis_0.9.5.0", "SanntiS 0.9.5.0", "*.sanntis.gff", sanntis09_regions),
         ("bgc-prophet_0.1.2", "BGC-Prophet 0.1.2", "prophet.csv", None)]  # parser needs the CDS tables


def intervals(tool_dirs: list[Path], name: str, cds: dict) -> dict[str, dict[str, list]]:
    """Tool run → record id → nucleotide intervals, for the runs present under the ``tool_dirs`` (a run found in a
    later directory replaces the same run of an earlier one: regenerated outputs over the archived ones)."""
    out: dict[str, dict[str, list]] = {}

    def collect(run: str, pattern: str, parse, tracks=TRACKS) -> None:
        for tools in tool_dirs:
            if not (tools / run).is_dir():
                continue
            files = [f for track in tracks for f in sorted((tools / run / f"{name}_{track}").glob(pattern))]
            if not files:  # an unfinished run (empty directories) is not a run without predictions
                log.warning("%s: no %s files under %s, left out", run, pattern, tools / run)
                continue
            d: dict[str, list] = {}
            for f in files:
                d.update(parse(f))
            out[run] = d

    for lv in LEVELS:
        collect(f"antismash_7.1.0_{lv}", "*.json.gz", T.antismash_regions)
    collect("plantismash_1.0", "*.gbk.gz", T.plantismash_regions, tracks=("euk",))
    tables = {cid: c for cid, (_, c) in cds.items()}
    for run, _, pattern, parse in PLAIN:
        collect(run, pattern, parse or (lambda f: prophet_regions(f, tables)))
    return out


def scopes(s: M.Scored, ids: set[str], *, shuffled: bool, n_boot: int) -> dict:
    """Subset metrics with bootstrap CIs (over the subset's contigs): ``proka`` / ``euk`` (every prediction and
    BGC on contigs of that kingdom) and, on the shuffled sets only, ``fungi``: the contigs whose one test BGC is a
    fungal BGC of the pool, so precision and recall are both about fungal BGCs. The isolate contigs carry several
    BGCs of mixed lineages, so their predictions cannot be attributed to a lineage: no ``fungi`` scope there (its
    fungal recall is in the strata / radar tables)."""
    res = {}
    for key, (sub, mask) in subsets(s, ids, shuffled=shuffled).items():
        res[f"metrics_{key}"] = {**M.aggregate(sub, gt_mask=mask),
                                 "n_contigs": float(pd.concat([sub.gt["unit"], sub.pred["unit"]]).nunique())}
        res[f"ci_{key}"] = SM.bootstrap(sub, gt_mask=mask, n_boot=n_boot)
    return res


def subsets(s: M.Scored, ids: set[str], *, shuffled: bool) -> dict[str, tuple[M.Scored, np.ndarray]]:
    """The scopes of ``scopes()``: scope → (subset, ground-truth mask); empty scopes left out."""
    out = {}
    for k in ("proka", "euk"):
        sub = M.Scored(s.gt[s.gt["kingdom"] == k], s.pred[s.pred["kingdom"] == k])
        out[k] = (sub, sub.gt["bgc"].isin(ids).to_numpy())
    if shuffled:
        if (n := s.gt.groupby("unit")["bgc"].nunique().max()) > 1:
            raise ValueError(f"shuffled contig with {n} test BGCs: the fungal scope needs one BGC per contig")
        units = set(s.gt.loc[s.gt["bgc"].isin(ids) & s.gt["lineage"].eq("fungi"), "unit"])
        sub = M.Scored(s.gt[s.gt["unit"].isin(units)], s.pred[s.pred["unit"].isin(units)])
        out["fungi"] = (sub, np.ones(len(sub.gt), bool))
    return {k: v for k, v in out.items() if v[1].any()}


PAIRED = ("SanntiS", "antiSMASH 7.1.0 relaxed (lineage-matched)")  # (a, b): differences are a − b
PAIRED_METRICS = ("precision_t", "recall_t", "f1_macro_t", "bgc_precision", "bgc_recall")


def paired(a: tuple[M.Scored, np.ndarray], b: tuple[M.Scored, np.ndarray], *, n_boot: int,
           seed: int = 0) -> list[dict]:
    """Paired bootstrap of a − b on one scope: each replicate resamples the scope's contigs once and scores both
    methods on that same resample. 95% percentile interval and a two-sided bootstrap p (twice the smaller tail share
    of replicates at or beyond 0)."""
    (sa, ma), (sb, mb) = a, b
    units = pd.unique(pd.concat([sa.gt["unit"], sa.pred["unit"], sb.gt["unit"], sb.pred["unit"]], ignore_index=True))
    code = {u: i for i, u in enumerate(units)}
    sides = []
    for s, m in ((sa, ma), (sb, mb)):
        gt = s.gt[m]
        sides.append((SM._arrays(gt, s.pred), gt["unit"].map(code).to_numpy(int),
                      s.pred["unit"].map(code).to_numpy(int)))
    rng = np.random.default_rng(seed)
    d: dict[str, list[float]] = {x: [] for x in PAIRED_METRICS}
    for _ in range(n_boot):
        w = np.bincount(rng.integers(0, len(units), len(units)), minlength=len(units)).astype(float)
        va, vb = (SM._metrics(*arr, w[gc], w[pc]) for arr, gc, pc in sides)
        for x in PAIRED_METRICS:
            d[x].append(va[x] - vb[x])
    pa, pb = M.aggregate(sa, gt_mask=ma), M.aggregate(sb, gt_mask=mb)
    rows = []
    for x in PAIRED_METRICS:
        v = np.asarray(d[x], float)
        v = v[~np.isnan(v)]
        rows.append({"metric": x, "a": pa[x], "b": pb[x], "diff": pa[x] - pb[x],
                     "ci_lo": float(np.percentile(v, 2.5)), "ci_hi": float(np.percentile(v, 97.5)),
                     "p": float(min(1.0, 2 * min((v <= 0).mean(), (v >= 0).mean()))), "n_boot": len(v),
                     "n_contigs": len(units)})
    return rows


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--testsets", type=Path, required=True)
    p.add_argument("--tools", type=Path, action="append", default=[])
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--test-bgcs", type=Path, required=True)
    p.add_argument("--tokens", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--n-boot", type=int, default=1000)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    a.out.mkdir(parents=True, exist_ok=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = SanntisModel.from_pretrained(a.model)
    t_star = model.config.decode.threshold
    tb = pd.read_parquet(a.test_bgcs)
    meta = tb.set_index("entity_id").assign(classes=lambda d: d["classes"].map(tuple))
    pools = {"A": set(tb.loc[tb["test_a"], "entity_id"]), "B": set(tb.loc[tb["test_b"], "entity_id"])}
    plants = set(meta.index[meta["lineage"] == "plants"])
    tok = pd.read_parquet(a.tokens)
    tokens = dict(zip(tok["md5"], map(list, tok["tokens"]), strict=True))
    grid = sorted(set(M.threshold_grid(load_recipe().operating_point.grid)) | {t_star})

    results, rows, prows = {}, [], []
    for name in SETS:
        kind = KIND[name]
        cds = load_cds(a.testsets, name)
        cs = contigs(cds, tokens)
        orf = {cid: (c["start"].to_numpy(), c["end"].to_numpy(),
                     [(int(g["start"].min()), int(g["end"].max())) for e, g in c.groupby("entity_id") if e in plants]
                     if track == "euk" else []) for cid, (track, c) in cds.items()}
        ivs = intervals(a.tools, name, cds)
        regions = {}
        for lv in LEVELS:
            run = f"antismash_7.1.0_{lv}"
            if run in ivs:
                regions[f"antiSMASH 7.1.0 {lv} (lineage-matched)"] = {
                    cid: T.to_orf_regions(s, e, T.lineage_matched(ivs[run].get(cid, []),
                                                                  ivs.get("plantismash_1.0", {}).get(cid, []), spans))
                    for cid, (s, e, spans) in orf.items()}
        for run, label, _, _ in PLAIN:
            if run in ivs:
                regions[label] = {cid: T.to_orf_regions(s, e, ivs[run].get(cid, [])) for cid, (s, e, _) in orf.items()}
        preds = predict_contigs(model, cs, grid)
        for test in ("A", "B"):
            key = f"{test}-{kind}"
            ids = pools[test]
            ours = evaluate_dataset(preds, t_star, grid, meta, ids, n_boot=a.n_boot, shuffled=kind == "shuffled")
            sub = {"shuffled": kind == "shuffled", "n_boot": a.n_boot}
            sc = {"SanntiS": scored(preds, t_star, meta)}
            res = {"SanntiS": {**ours, **scopes(sc["SanntiS"], ids, **sub)}}
            for tool, reg in regions.items():
                tp = T.tool_preds(cs, reg)
                r = {k: v for k, v in T.score_tool(tp, meta, ids, n_boot=a.n_boot).items()
                     if not k.startswith("metrics_")}  # replaced by scopes(): CIs, fungi on shuffled sets only
                sc[tool] = scored(tp, 1.0, meta)
                res[tool] = {**r, **scopes(sc[tool], ids, **sub)}
            if all(t in sc for t in PAIRED):
                sa, sb = (subsets(sc[t], ids, shuffled=kind == "shuffled") for t in PAIRED)
                for scope in ("euk", "fungi"):
                    if scope in sa:
                        prows += [{"set": key, "scope": scope, "a": PAIRED[0], "b": PAIRED[1], **r}
                                  for r in paired(sa[scope], sb[scope], n_boot=a.n_boot)]
            results[key] = res
            for tool, r in res.items():
                rows.append({"set": key, "tool": tool, "n_bgc": int(r["metrics"]["n_gt"]),
                             **{k: round(r["metrics"][k], 4) for k in M.METRICS},
                             **{f"f1_macro_t_{k}": round(r[f"metrics_{k}"]["f1_macro_t"], 4)
                                for k in ("proka", "euk", "fungi") if f"metrics_{k}" in r}})
                log.info("%s %s F1_macro_T %.3f", key, tool, r["metrics"]["f1_macro_t"])
    order = ["A-isolate", "A-shuffled", "B-isolate", "B-shuffled"]
    results = {k: results[k] for k in order}
    (a.out / "results.pkl").write_bytes(pickle.dumps({"t_star": t_star, "results": results}))
    pd.DataFrame(rows).to_csv(a.out / "summary.tsv", sep="\t", index=False)
    if prows:  # SanntiS 1.0 − antiSMASH relaxed, paired over contigs, on the eukaryotic / fungal scopes
        pd.DataFrame(prows).to_csv(a.out / "paired.tsv", sep="\t", index=False, float_format="%.6g")
    expected.write(results, a.out)
    log.info("wrote %s", a.out)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
