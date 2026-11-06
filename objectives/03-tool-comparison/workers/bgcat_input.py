"""BGCat input as its README prescribes for antiSMASH users: per test BGC, the proteins of the antiSMASH genes labelled
``biosynthetic`` / ``biosynthetic-additional`` inside an antiSMASH region and inside the BGC.

    python bgcat_input.py --testsets DIR --antismash DIR --out bgcs.json

``--antismash``: an antiSMASH 7.1.0 run on the nucleotide test sets (``<set>_<track>/<contig>.json.gz``); the isolate
set holds every test BGC once. BGCs without such genes get no BGCat call (reported as coverage).
"""

import argparse
import gzip
import json
import re
from pathlib import Path

from testsets import TRACKS, load_cds

KINDS = {"biosynthetic", "biosynthetic-additional"}
SET = "isolate"


def span(loc: str) -> tuple[int, int]:
    """antiSMASH JSON location string (``[a:b](+)``, ``join{[a:b](+), ...}``) → outer 0-based span."""
    xs = [int(v) for pair in re.findall(r"\[<?>?(\d+):<?>?(\d+)\]", loc) for v in pair]
    return min(xs), max(xs)


p = argparse.ArgumentParser()
p.add_argument("--testsets", type=Path, required=True)
p.add_argument("--antismash", type=Path, required=True)
p.add_argument("--out", type=Path, required=True)
a = p.parse_args()

cds = load_cds(a.testsets, SET)
bgcs: dict[str, list[str]] = {}
for track in TRACKS:
    for f in sorted((a.antismash / f"{SET}_{track}").glob("*.json.gz")):
        with gzip.open(f, "rt") as fh:
            records = json.load(fh)["records"]
        for rec in records:
            _, c = cds[rec["id"]]
            pos = c[c["source_kind"] == "positive"]
            spans = {e: (int(g["start"].min()), int(g["end"].max())) for e, g in pos.groupby("entity_id")}
            areas = [(int(x["start"]), int(x["end"])) for x in rec.get("areas", [])]
            for feat in rec["features"]:
                q = feat["qualifiers"]
                if feat["type"] != "CDS" or not KINDS & set(q.get("gene_kind", [])) or "translation" not in q:
                    continue
                s, e = span(feat["location"])
                if not any(s < b and e > a_ for a_, b in areas):
                    continue
                for bgc, (bs, be) in spans.items():
                    if s >= bs and e <= be:
                        bgcs.setdefault(bgc, []).append(q["translation"][0])
a.out.write_text(json.dumps(bgcs))
n_bgc = sum(len(set(c.loc[c["source_kind"] == "positive", "entity_id"])) for _, c in cds.values())
print(f"{len(bgcs)} of {n_bgc} test BGCs with antiSMASH biosynthetic genes -> {a.out}")
