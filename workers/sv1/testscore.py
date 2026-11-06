# Copied from SanntiS 842675a `src/sanntis/train/evaluate.py` + `TestContig` / `select_test_pools` of `testsets.py`
# Test-set code leaves SanntiS in EMG-10233; only the imports were changed.
"""Test-set scoring for ``--evaluate`` (spec §7, §8.3): a bundle on the isolate / shuffled test contigs.

One detector pass per contig; decode + RF filter at every sweep threshold; metrics at t* with bootstrap CIs over
source contigs; the four PR curves and their AUPRC; strata (class, lineage, novelty) as GT masks. Test A is
scored on the Test B contigs: its recall reads the active BGCs only, precision is shared (pending BGCs stay
ground truth). Shuffled sets also report the native-BGC sensitivity row (natives counted as background).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sanntis.features.tokens import linearize
from sanntis.model.bundle import SanntisModel
from sanntis.model.decode import decode, orf_scores
from sanntis.model.detector import predict_tokens

from sv1 import metrics as M
from sv1.novelty import BINS

V3_MAX = "BGC0002750"


@dataclass
class TestContig:
    __test__ = False  # not a pytest class

    contig_id: str
    kingdom: str
    cds_tokens: list[list[str]]
    orf_bgc: list[str | None]
    native: np.ndarray
    source: dict = field(default_factory=dict)


def select_test_pools(positives: pd.DataFrame, *, max_locus_bp: int = 200_000) -> dict[str, set[str]]:
    """``positives``: entity_id, status (5.0rc1), length, cds_tokens → {"A": ids, "B": ids, "excluded": ids}."""
    p = positives[positives["entity_id"] >= V3_MAX]
    bad = p[(p["length"] > max_locus_bp) | (p["cds_tokens"].map(len) == 0)]
    p = p[~p["entity_id"].isin(bad["entity_id"])]
    return {"A": set(p.loc[p["status"] == "active", "entity_id"]),
            "B": set(p.loc[p["status"].isin(["active", "pending"]), "entity_id"]),
            "excluded": set(bad["entity_id"])}


log = logging.getLogger(__name__)


@dataclass
class ContigPred:
    contig: TestContig
    scores: np.ndarray
    regions: dict[float, list[tuple[int, int, str, float]]]  # threshold → kept (a, b, class, p_non)


def predict_contigs(model: SanntisModel, contigs: list[TestContig], thresholds: list[float]) -> list[ContigPred]:
    f = model.config.features
    dets = model.detectors("cpu")
    w = model.config.detector.window
    out = []
    for c in contigs:
        ts = linearize(c.cds_tokens, model.vocab, unk_policy=f.unk_policy, gene_sep=f.gene_sep,
                       dedup_contiguous=f.dedup_contiguous)
        p = np.mean([predict_tokens(m, ts.tokens, w) for m in dets], axis=0)
        s = orf_scores(p, ts.orf_index, len(c.cds_tokens), agg=model.config.decode.agg,
                       fill=model.config.decode.fill)
        reg = {}
        for t in thresholds:
            r = decode(s, model.config.decode, threshold=t)
            kept = []
            if r:
                keep, cls, pn = model.postfilter.classify([c.cds_tokens[a:b] for a, b in r],
                                                          model.config.postfilter_drop_threshold)
                kept = [(a, b, cl, float(q)) for (a, b), k, cl, q in zip(r, keep, cls, pn, strict=True) if k]
            reg[t] = kept
        out.append(ContigPred(c, s, reg))
    return out


def scored(preds: list[ContigPred], t: float, bgc_meta: pd.DataFrame, *, natives: bool = True) -> M.Scored:
    """``bgc_meta`` indexed by entity_id: classes (tuple), lineage, novelty, kingdom."""
    parts = []
    for cp in preds:
        c = cp.contig
        orf_bgc = np.array(c.orf_bgc, dtype=object)
        classes = {b: tuple(bgc_meta.at[b, "classes"]) for b in set(x for x in c.orf_bgc if x is not None)}
        s = M.score_contig([(a, b) for a, b, *_ in cp.regions[t]], orf_bgc, unit=c.contig_id,
                           orf_native=c.native if natives and c.native.any() else None, bgc_classes=classes)
        s.gt["kingdom"] = c.kingdom
        s.pred["kingdom"] = c.kingdom
        parts.append(s)
    s = M.concat(parts)
    for col in ("lineage", "novelty"):
        s.gt[col] = s.gt["bgc"].map(bgc_meta[col]) if col in bgc_meta else None
    return s


LINEAGES = ("prokaryotes", "fungi", "plants", "other eukaryotes")


def radar_table(s: M.Scored, gt_mask: np.ndarray | None = None) -> pd.DataFrame:
    """Recall_T and BGC recall per class, lineage and novelty stratum (one row per axis; empty strata left out)."""
    mask = np.ones(len(s.gt), bool) if gt_mask is None else np.asarray(gt_mask, bool)
    a = M.aggregate(s, gt_mask=mask)
    rows = [{"sector": "class", "axis": c, "n": int(round(a[f"n_gt[{c}]"])), "recall_t": a[f"recall_t[{c}]"],
             "bgc_recall": a[f"bgc_recall[{c}]"]} for c in M.CLASSES if a.get(f"n_gt[{c}]", 0) > 0]
    for col, values in (("lineage", LINEAGES), ("novelty", BINS)):
        for v in values if col in s.gt else ():
            m = mask & (s.gt[col] == v).to_numpy()
            if m.any():
                b = M.aggregate(s, gt_mask=m)
                rows.append({"sector": col, "axis": v, "n": int(m.sum()), "recall_t": b["recall_t"],
                             "bgc_recall": b["bgc_recall"]})
    return pd.DataFrame(rows, columns=["sector", "axis", "n", "recall_t", "bgc_recall"])


def ignored_predictions(preds: list[ContigPred], t: float) -> int:
    n = 0
    for cp in preds:
        lab = np.array([x is not None for x in cp.contig.orf_bgc])
        for a, b, *_ in cp.regions[t]:
            n += bool(cp.contig.native[a:b].any() and not lab[a:b].any())
    return n


def pr_curves(preds: list[ContigPred], thresholds: list[float], bgc_meta: pd.DataFrame,
              gt_ids: set[str]) -> pd.DataFrame:
    rows = []
    for t in thresholds:
        s = scored(preds, t, bgc_meta)
        a = M.aggregate(s, gt_mask=s.gt["bgc"].isin(gt_ids).to_numpy())
        rows.append({"threshold": t, **{k: a[k] for k in M.METRICS}})
    return pd.DataFrame(rows)


CURVES = {"PrecT_vs_RecMacroT": ("recall_macro_t", "precision_t"), "PrecT_vs_RecT": ("recall_t", "precision_t"),
          "BGCprec_vs_BGCrec": ("bgc_recall", "bgc_precision"),
          "BGCprec_vs_BGCmacrorec": ("bgc_macro_recall", "bgc_precision")}


def auprcs(curve: pd.DataFrame) -> dict[str, float]:
    return {name: M.auprc(curve[x].to_numpy(), curve[y].to_numpy()) for name, (x, y) in CURVES.items()}


def evaluate_dataset(preds: list[ContigPred], t_star: float, thresholds: list[float], bgc_meta: pd.DataFrame,
                     gt_ids: set[str], *, n_boot: int = 1000, shuffled: bool = False) -> dict:
    s = scored(preds, t_star, bgc_meta)
    mask = s.gt["bgc"].isin(gt_ids).to_numpy()
    res = {"metrics": M.aggregate(s, gt_mask=mask), "ci": M.bootstrap(s, gt_mask=mask, n_boot=n_boot)}
    for k in ("proka", "euk"):
        sub = M.Scored(s.gt[s.gt["kingdom"] == k], s.pred[s.pred["kingdom"] == k])
        m = sub.gt["bgc"].isin(gt_ids).to_numpy()
        if m.any():
            res[f"metrics_{k}"] = M.aggregate(sub, gt_mask=m)
    strata = []
    for col in ("lineage", "novelty"):
        for v in sorted(x for x in s.gt[col].dropna().unique()) if col in s.gt else []:
            m = mask & (s.gt[col] == v).to_numpy()
            if m.any():
                a = M.aggregate(s, gt_mask=m)
                strata.append({"axis": col, "stratum": v, "n": int(m.sum()), "bgc_macro_recall": a["bgc_macro_recall"],
                               "recall_macro_t": a["recall_macro_t"]})
    res["strata"] = pd.DataFrame(strata)
    res["radar"] = radar_table(s, mask)
    curve = pr_curves(preds, thresholds, bgc_meta, gt_ids)
    res["curve"] = curve
    res["auprc"] = auprcs(curve)
    if shuffled:
        res["ignored_native_predictions"] = ignored_predictions(preds, t_star)
        s2 = scored(preds, t_star, bgc_meta, natives=False)
        a2 = M.aggregate(s2, gt_mask=s2.gt["bgc"].isin(gt_ids).to_numpy())
        res["natives_as_background"] = {"precision_t": a2["precision_t"], "bgc_precision": a2["bgc_precision"]}
    return res


# ── NP classification on the test sets (spec §7.5, §8.4) ────────────────────


def npc_scores(pp: np.ndarray, ps: np.ndarray, yp: np.ndarray, ys: np.ndarray, m: np.ndarray, pathways: list[str],
               abstention: dict, targets=(0.9, 0.8)) -> tuple[dict, pd.DataFrame]:
    """``pp`` remapped pathway probabilities (calls at >= 0.5, or the top one), ``ps`` superclass probabilities."""
    from sklearn.metrics import average_precision_score

    calls = pp >= 0.5
    calls[np.arange(len(pp)), pp.argmax(1)] |= ~calls.any(1)
    per = []
    for j, name in enumerate(pathways):
        y = yp[:, j] > 0.5
        if not y.any():
            continue
        f1 = 2 * (y & calls[:, j]).sum() / max(y.sum() + calls[:, j].sum(), 1)
        per.append({"pathway": name, "ap": float(average_precision_score(y, pp[:, j])), "f1": float(f1),
                    "support": int(y.sum())})
    per = pd.DataFrame(per, columns=["pathway", "ap", "f1", "support"])
    top1 = yp[np.arange(len(pp)), pp.argmax(1)] > 0.5
    sm = m > 0.5
    out = {"n": len(pp), "pathway_macro_ap": float(per["ap"].mean()) if len(per) else np.nan,
           "pathway_macro_f1": float(per["f1"].mean()) if len(per) else np.nan,
           "pathway_top1": float(top1.mean()) if len(pp) else np.nan, "n_superclass": int(sm.sum())}
    if sm.any():
        Y, P = ys[sm], ps[sm]
        order = np.argsort(-P, 1)
        out["superclass_top1"] = float((Y[np.arange(len(Y)), order[:, 0]] > 0.5).mean())
        out["superclass_top3"] = float(np.mean([Y[i, order[i, :3]].max() > 0.5 for i in range(len(Y))]))
        conf, correct = P.max(1), Y[np.arange(len(Y)), order[:, 0]] > 0.5
        for x in targets:
            ans = conf >= abstention["sc"][x]
            out[f"coverage@{x}"] = float(ans.mean())
            out[f"precision@{x}"] = float(correct[ans].mean()) if ans.any() else np.nan
    return out, per


def npc_evaluate(model, preds: list[ContigPred], t: float, tests: pd.DataFrame, labels: pd.DataFrame,
                 parents: dict, gt_ids: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """True loci (ceiling) and predicted regions at ``t`` of the BGCs in ``gt_ids`` with an NPC label.
    A region takes the labels of the BGC with which it shares the most ORFs."""
    from sanntis.train.npc_train import LabelSpace, labelled

    npc = model.npc
    ls = LabelSpace(superclasses=list(npc.superclasses), parents=parents, pathways=list(npc.pathways))
    ok = set(labelled(labels, sorted(gt_ids)))
    rows, pers = [], []
    t_idx = tests.set_index("entity_id")
    ids = sorted(ok & set(t_idx.index))
    loci = {"cds": [list(t_idx.at[i, "cds_tokens"]) for i in ids], "kingdom": [t_idx.at[i, "kingdom"] for i in ids],
            "ids": ids}
    reg = {"cds": [], "kingdom": [], "ids": []}
    for cp in preds:
        c = cp.contig
        bgc = np.array([x if x is not None else "" for x in c.orf_bgc], dtype=object)
        for a, b, *_ in cp.regions[t]:
            ents, cnt = np.unique(bgc[a:b][bgc[a:b] != ""], return_counts=True)
            if not len(ents):
                continue
            main = str(ents[int(np.argmax(cnt))])
            if main in ok:
                reg["cds"].append(c.cds_tokens[a:b])
                reg["kingdom"].append(c.kingdom)
                reg["ids"].append(main)
    for name, ex in (("true loci", loci), ("predicted regions", reg)):
        if not ex["ids"]:
            continue
        pp, ps = npc.proba(ex["cds"], np.asarray(ex["kingdom"]))
        yp, ys, m = ls.encode(labels, ex["ids"])
        s, per = npc_scores(pp, ps, yp, ys, m, list(npc.pathways), npc.thresholds.abstention)
        rows.append({"input": name, "stratum": "all", **s})
        fungi = np.array([str(t_idx.at[i, "lineage"]) == "fungi" if "lineage" in t_idx else False for i in ex["ids"]])
        if fungi.any():
            sf, _ = npc_scores(pp[fungi], ps[fungi], yp[fungi], ys[fungi], m[fungi], list(npc.pathways),
                               npc.thresholds.abstention)
            rows.append({"input": name, "stratum": "fungi", **sf})
        pers.append(per.assign(input=name))
    return pd.DataFrame(rows), (pd.concat(pers, ignore_index=True) if pers else pd.DataFrame())
