"""Paired bootstrap of SanntiS 1.2.0 − each tool on the BGC-level metrics of B-shuffled (D3.2 headline numbers).

    python paired_bgc.py --testsets DIR --tools DIR --model BUNDLE --test-bgcs PARQUET --tokens PARQUET --out TSV

Same inputs, regions and scoring as ``score.py`` (lineage-matched antiSMASH, every tool on the same contigs and CDS),
restricted to Test B on the shuffled set. Scopes: ``all`` (every contig), ``proka`` / ``euk`` (contigs of that
kingdom) and ``fungi`` (contigs whose one test BGC is fungal), as in ``score.subsets``. Each replicate resamples the
scope's contigs once and scores both methods on that resample (``score.paired``).
"""
import argparse
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import score as S
from sanntis.model.bundle import SanntisModel
from sanntis.train import metrics as M
from testsets import KIND, contigs, load_cds

from sv1 import tools as T
from sv1.testscore import predict_contigs, scored

METRICS = ("f1_macro", "bgc_precision", "bgc_recall", "bgc_macro_recall", "f1_macro_t")
SET = "shuffled"
log = logging.getLogger("paired_bgc")


def main() -> None:
    p = argparse.ArgumentParser()
    for f in ("--testsets", "--model", "--test-bgcs", "--tokens", "--out"):
        p.add_argument(f, type=Path, required=True)
    p.add_argument("--tools", type=Path, action="append", default=[])
    p.add_argument("--n-boot", type=int, default=1000)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    S.PAIRED_METRICS = METRICS  # score.paired reads the module constant

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = SanntisModel.from_pretrained(a.model)
    t_star = model.config.decode.threshold
    tb = pd.read_parquet(a.test_bgcs)
    meta = tb.set_index("entity_id").assign(classes=lambda d: d["classes"].map(tuple))
    ids = set(tb.loc[tb["test_b"], "entity_id"])
    plants = set(meta.index[meta["lineage"] == "plants"])
    tok = pd.read_parquet(a.tokens)
    tokens = dict(zip(tok["md5"], map(list, tok["tokens"]), strict=True))

    name = next(n for n in KIND if KIND[n] == SET)
    cds = load_cds(a.testsets, name)
    cs = contigs(cds, tokens)
    orf = {cid: (c["start"].to_numpy(), c["end"].to_numpy(),
                 [(int(g["start"].min()), int(g["end"].max())) for e, g in c.groupby("entity_id") if e in plants]
                 if track == "euk" else []) for cid, (track, c) in cds.items()}
    ivs = S.intervals(a.tools, name, cds)
    regions = {}
    run = "antismash_7.1.0_relaxed"
    if run in ivs:
        regions["antiSMASH 7.1.0 relaxed (lineage-matched)"] = {
            cid: T.to_orf_regions(s, e, T.lineage_matched(ivs[run].get(cid, []),
                                                          ivs.get("plantismash_1.0", {}).get(cid, []), spans))
            for cid, (s, e, spans) in orf.items()}
    for run, label, _, _ in S.PLAIN:
        if run in ivs and label != "GECCO 0.9.5":
            regions[label] = {cid: T.to_orf_regions(s, e, ivs[run].get(cid, [])) for cid, (s, e, _) in orf.items()}

    preds = predict_contigs(model, cs, [t_star])
    sc = {"SanntiS": scored(preds, t_star, meta)}
    for tool, reg in regions.items():
        sc[tool] = scored(T.tool_preds(cs, reg), 1.0, meta)

    def scopes(s: M.Scored) -> dict:
        out = {"all": (s, s.gt["bgc"].isin(ids).to_numpy())}
        out.update(S.subsets(s, ids, shuffled=True))
        return out

    ours = scopes(sc["SanntiS"])
    rows = []
    for tool in regions:
        theirs = scopes(sc[tool])
        for scope in ("all", "proka", "euk", "fungi"):
            if scope in ours and scope in theirs:
                rows += [{"set": f"B-{SET}", "scope": scope, "method_a": "SanntiS 1.2.0", "method_b": tool, **r}
                         for r in S.paired(ours[scope], theirs[scope], n_boot=a.n_boot)]
                log.info("%s %s done", tool, scope)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False, float_format="%.6g")
    log.info("wrote %s", a.out)


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
