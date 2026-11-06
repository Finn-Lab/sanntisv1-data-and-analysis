"""Protein FASTA of one test contig in gene order. Headers ``<record>__<cds index>`` (BGC-Prophet input), or with
``cds_id`` the CDS id (= the GenBank ``/protein_id`` / ``/locus_tag``; InterProScan input for SanntiS 0.9.5).

    python proteins.py <contig.cds.parquet> <out.faa> [cds_id]
"""

import sys
from pathlib import Path

import pandas as pd

c = pd.read_parquet(sys.argv[1]).drop_duplicates(["start", "end"]).sort_values("start").reset_index(drop=True)
by_id = len(sys.argv) > 3 and sys.argv[3] == "cds_id"
with open(Path(sys.argv[2]), "w") as f:
    for i, (cid, pid, s) in enumerate(zip(c["contig_id"], c["cds_id"], c["translation"], strict=True)):
        f.write(f">{pid if by_id else f'{cid}__{i}'}\n{s.rstrip('*').replace('J', 'X')}\n")
