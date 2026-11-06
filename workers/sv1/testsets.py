"""The nucleotide test sets on disk (objective 02's ``testsets/`` tree, the ``nucleotide_testsets`` deposits):
``<root>/<set>/<track>/<prefix>.{gbk,cds.parquet,ground_truth.bed}`` → CDS tables and test contigs."""

from pathlib import Path

import pandas as pd
from sanntis.io.fasta import clean_protein, md5

from sv1.testscore import TestContig

TRACKS = ("proka", "euk")


def sets_of(root: Path) -> list[str]:
    """The set directories of a test-set tree (sorted)."""
    return sorted(d.name for d in Path(root).iterdir() if d.is_dir() and any((d / t).is_dir() for t in TRACKS))


def load_cds(root: Path, name: str) -> dict[str, tuple[str, pd.DataFrame]]:
    """Record id → (track, CDS table sorted by start, one row per distinct span, with ``md5`` of the protein)."""
    out = {}
    for track in TRACKS:
        for f in sorted((Path(root) / name / track).glob("*.cds.parquet")):
            c = pd.read_parquet(f).drop_duplicates(["start", "end"]).sort_values("start").reset_index(drop=True)
            if "native_bgc" not in c:
                c["native_bgc"] = False
            c["md5"] = [md5(clean_protein(s)) for s in c["translation"]]
            out[c["contig_id"].iloc[0]] = (track, c)
    return out


def proteins_of(root: Path, names: list[str] | None = None) -> dict[str, str]:
    """md5 → cleaned protein of every CDS of the given sets (default: all)."""
    return {m: clean_protein(s) for name in (names or sets_of(root)) for _, c in load_cds(root, name).values()
            for m, s in zip(c["md5"], c["translation"], strict=True)}


def contigs(cds: dict[str, tuple[str, pd.DataFrame]], tokens: dict[str, list[str]]) -> list[TestContig]:
    """Test contigs with SanntiS tokens per CDS (``tokens``: protein md5 → token list)."""
    out = []
    for cid, (track, c) in cds.items():
        bgc = [e if k == "positive" else None for k, e in zip(c["source_kind"], c["entity_id"], strict=True)]
        out.append(TestContig(cid, track, [tokens.get(m, []) for m in c["md5"]], bgc,
                              c["native_bgc"].fillna(False).astype(bool).to_numpy()))
    return out
