"""SanntiS tokens and Pfam sets of every test-contig protein: InterProScan 6 through SanntiS's annotation (its
md5-keyed cache, so proteins already annotated for training are not run again).

    python tokens.py --testsets DIR --out tokens.parquet [--container-engine E] [--nf-config F ...]

Writes one row per protein md5: ``tokens`` (the InterPro signature map of the SanntiS recipe, in hit order) and
``pfam`` (the Pfam accessions hit, for the novelty strata). The IPS6 run directory sits next to the output.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
from sanntis.annotate.ips6 import annotate, default_datadir
from sanntis.config import NextflowOptions, load_recipe
from sanntis.io.ips_tsv import signature_map_from_entries, tokens_by_key

from sv1.testsets import proteins_of

p = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
p.add_argument("--testsets", type=Path, required=True)
p.add_argument("--out", type=Path, required=True)
p.add_argument("--container-engine", default=None)
p.add_argument("--nf-config", type=Path, action="append", default=[])
a = p.parse_args()
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

proteins = proteins_of(a.testsets)
logging.info("%d unique proteins", len(proteins))
ann = annotate(proteins, NextflowOptions(container_engine=a.container_engine, nf_config=a.nf_config),
               a.out.parent / "ips6_run")
logging.info("annotated: %d cached, %d sent to InterProScan 6", ann.n_cached, ann.n_sent)
sm = signature_map_from_entries(default_datadir() / "interpro" / load_recipe().ips.interpro / "entries.json")
tok = tokens_by_key(ann.hits, sm, key="md5")
pfam = ann.hits[ann.hits["app"] == "pfam"].groupby("md5")["signature"].agg(lambda s: sorted(set(s))).to_dict()
a.out.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame({"md5": list(proteins), "tokens": [list(tok.get(m, [])) for m in proteins],
              "pfam": [list(pfam.get(m, [])) for m in proteins]}).to_parquet(a.out, index=False)
logging.info("%d / %d proteins with tokens -> %s", sum(bool(tok.get(m)) for m in proteins), len(proteins), a.out)
