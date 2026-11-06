# Copied from SanntiS 842675a `src/sanntis/train/baselines.py`
# Test-set code leaves SanntiS in EMG-10233; only the imports were changed.
"""Baseline tools on the test sets (spec §8.3): antiSMASH, plantiSMASH and GECCO predictions scored with the same
ORF-space metrics as SanntiS.

The tools need nucleotide contigs, so baselines are scored on nucleotide-built test sets whose CDS table gives each
CDS's coordinates; their outputs are parsed here (no tool dependency) and projected to ORF regions. antiSMASH is
lineage-matched as a user would run it: ``--taxon bacteria`` on prokaryotic contigs, ``--taxon fungi`` (fungiSMASH)
on eukaryotic ones, and plantiSMASH for regions over a plant BGC (plant BGCs sit on fungal / protist backbones, so the
mode cannot be chosen per contig).
"""

from __future__ import annotations

import gzip
import json
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from sv1 import metrics as M
from sv1.testscore import ContigPred, radar_table, scored

Interval = tuple[int, int, float]  # 0-based half-open nucleotide span, score


def _open(path: Path):
    return gzip.open(path, "rt") if path.name.endswith(".gz") else open(path)


def antismash_regions(path: Path) -> dict[str, list[Interval]]:
    """antiSMASH JSON: ``records[].areas`` (regions) per record id."""
    with _open(path) as fh:
        data = json.load(fh)
    return {r["id"]: [(int(a["start"]), int(a["end"]), 1.0) for a in r.get("areas", [])]
            for r in data.get("records", [])}


def plantismash_regions(path: Path) -> dict[str, list[Interval]]:
    """plantiSMASH 1.0 GenBank ``cluster`` features; record ids restored through ``<prefix>.ids.tsv`` if present."""
    from Bio import SeqIO

    ids_tsv = path.with_name(path.name.split(".gbk")[0] + ".ids.tsv")
    ids = dict(line.split("\t") for line in ids_tsv.read_text().splitlines() if line) if ids_tsv.exists() else {}
    out: dict[str, list[Interval]] = {}
    with _open(path) as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            rid = ids.get(rec.id, ids.get(rec.name, rec.id))
            out.setdefault(rid, []).extend((int(f.location.start), int(f.location.end), 1.0)
                                           for f in rec.features if f.type == "cluster")
    return out


def gecco_regions(path: Path) -> dict[str, list[Interval]]:
    """GECCO ``*.clusters.tsv`` (1-based closed coordinates), scored by ``average_p``."""
    out: dict[str, list[Interval]] = {}
    for r in pd.read_csv(path, sep="\t").itertuples(index=False):
        out.setdefault(str(r.sequence_id), []).append((int(r.start) - 1, int(r.end), float(r.average_p)))
    return out


def to_orf_regions(starts: np.ndarray, ends: np.ndarray, ivs: Sequence[Interval]) -> list[tuple[int, int]]:
    """Nucleotide intervals → ``[first, last + 1)`` over the CDS they overlap; overlapping projections merge."""
    out = []
    for a, b, _ in sorted(ivs):
        idx = np.flatnonzero((starts < b) & (ends > a))
        if idx.size:
            out.append((int(idx[0]), int(idx[-1]) + 1))
    merged: list[tuple[int, int]] = []
    for r in out:
        if merged and r[0] < merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], r[1]))
        else:
            merged.append(r)
    return merged


def lineage_matched(base: Sequence[Interval], plant: Sequence[Interval],
                    plant_spans: Sequence[tuple[int, int]]) -> list[Interval]:
    """antiSMASH regions of one contig, with those over a plant BGC span replaced by the plantiSMASH regions there."""
    if not plant_spans:
        return list(base)
    hit = lambda iv: any(iv[0] < b and iv[1] > a for a, b in plant_spans)  # noqa: E731
    return [iv for iv in base if not hit(iv)] + [iv for iv in plant if hit(iv)]


def tool_preds(contigs: list, regions: dict[str, list[tuple[int, int]]]) -> list[ContigPred]:
    """Wrap one tool's ORF regions (contig id → regions) as single-threshold predictions (key 1.0)."""
    return [ContigPred(c, np.zeros(len(c.cds_tokens)),
                       {1.0: [(a, b, "", 0.0) for a, b in regions.get(c.contig_id, [])]}) for c in contigs]


def score_tool(preds: list[ContigPred], bgc_meta: pd.DataFrame, gt_ids: set[str], *, n_boot: int = 1000,
               t: float = 1.0) -> dict:
    """Metric set + CIs over contigs (as ``evaluate_dataset``), plus ``subset_metrics``."""
    s = scored(preds, t, bgc_meta)
    mask = s.gt["bgc"].isin(gt_ids).to_numpy()
    return {"metrics": M.aggregate(s, gt_mask=mask), "ci": M.bootstrap(s, gt_mask=mask, n_boot=n_boot),
            "radar": radar_table(s, mask), **subset_metrics(s, gt_ids)}


def subset_metrics(s: M.Scored, gt_ids: set[str]) -> dict:
    """``metrics_proka`` / ``metrics_euk`` (per kingdom) and ``metrics_fungi`` (recall over fungal BGCs, precision
    over every eukaryotic prediction)."""
    res = {}
    for k in ("proka", "euk"):
        sub = M.Scored(s.gt[s.gt["kingdom"] == k], s.pred[s.pred["kingdom"] == k])
        m = sub.gt["bgc"].isin(gt_ids).to_numpy()
        if m.any():
            res[f"metrics_{k}"] = M.aggregate(sub, gt_mask=m)
    if "lineage" in s.gt:
        sub = M.Scored(s.gt[s.gt["kingdom"] == "euk"], s.pred[s.pred["kingdom"] == "euk"])
        m = (sub.gt["bgc"].isin(gt_ids) & sub.gt["lineage"].eq("fungi")).to_numpy()
        if m.any():
            res["metrics_fungi"] = M.aggregate(sub, gt_mask=m)
    return res
