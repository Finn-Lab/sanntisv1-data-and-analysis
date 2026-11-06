"""Parsers of the tools added in this repository (the antiSMASH, plantiSMASH and GECCO parsers are in ``sv1.tools``).
Each returns record id → [(start, end, score)] nucleotide intervals (0-based, half-open)."""

import gzip
from pathlib import Path

import pandas as pd

Interval = tuple[int, int, float]


def _open(path: Path):
    return gzip.open(path, "rt") if path.name.endswith(".gz") else open(path)


def sanntis09_regions(path: Path) -> dict[str, list[Interval]]:
    """SanntiS 0.9 GFF3 ``CLUSTER`` features (1-based closed), scored by ``score=``."""
    out: dict[str, list[Interval]] = {}
    with _open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CLUSTER":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            out.setdefault(f[0], []).append((int(f[3]) - 1, int(f[4]), float(attrs.get("score", 1.0))))
    return out


def prophet_regions(path: Path, cds: dict[str, pd.DataFrame]) -> dict[str, list[Interval]]:
    """BGC-Prophet pipeline CSV: per 128-gene window, the BGC genes (``sentence``, headers ``<record>__<cds index>``)
    of windows called ``isBGC == Yes``; each becomes the span of those CDS (windows overlap; projections merge)."""
    out: dict[str, list[Interval]] = {}
    for r in pd.read_csv(path).itertuples(index=False):
        if r.isBGC != "Yes" or not isinstance(r.sentence, str):
            continue
        genes = [g.rsplit("__", 1) for g in r.sentence.split()]
        rec = genes[0][0]
        idx = [int(i) for _, i in genes]
        c = cds[rec]
        out.setdefault(rec, []).append((int(c["start"].iloc[min(idx)]), int(c["end"].iloc[max(idx)]), 1.0))
    return out
