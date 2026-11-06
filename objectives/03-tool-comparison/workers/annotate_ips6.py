"""SanntiS tokens of every test-set protein: InterProScan 6 through SanntiS's annotation (and its cache).

    python annotate_ips6.py --testsets DIR --out tokens.parquet [--container-engine E] [--nf-config F ...]

Writes one row per protein md5 with its token list (the InterPro 107.0 signature map of the SanntiS recipe).
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
from sanntis.annotate.ips6 import annotate, default_datadir
from sanntis.config import NextflowOptions, load_recipe
from sanntis.io.fasta import clean_protein
from sanntis.io.ips_tsv import signature_map_from_entries, tokens_by_key
from testsets import SETS, load_cds

p = argparse.ArgumentParser()
p.add_argument("--testsets", type=Path, required=True)
p.add_argument("--out", type=Path, required=True)
p.add_argument("--container-engine", default="docker")
p.add_argument("--nf-config", type=Path, action="append", default=[])
a = p.parse_args()
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

proteins = {m: clean_protein(s) for name in SETS for _, c in load_cds(a.testsets, name).values()
            for m, s in zip(c["md5"], c["translation"], strict=True)}
logging.info("%d unique proteins", len(proteins))
ann = annotate(proteins, NextflowOptions(container_engine=a.container_engine, nf_config=a.nf_config),
               a.out.parent / "ips6_run")
sm = signature_map_from_entries(default_datadir() / "interpro" / load_recipe().ips.interpro / "entries.json")
tok = tokens_by_key(ann.hits, sm, key="md5")
pd.DataFrame({"md5": list(proteins), "tokens": [list(tok.get(m, [])) for m in proteins]}).to_parquet(a.out)
logging.info("%d proteins with tokens -> %s", sum(bool(tok.get(m)) for m in proteins), a.out)
