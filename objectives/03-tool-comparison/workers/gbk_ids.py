"""antiSMASH / plantiSMASH need short record ids: shorten them in the input GenBank, restore them in the outputs.

    python gbk_ids.py shorten <in.gbk> <out.gbk> <ids.tsv>      # records renamed c1, c2, ...; ids.tsv: short → id
    python gbk_ids.py restore-json <antismash.json> <ids.tsv>   # in place: records[].id back to the original ids
"""

import json
import sys
from pathlib import Path

from Bio import SeqIO


def shorten(src: Path, dest: Path, ids: Path) -> None:
    recs, rows = [], []
    for i, r in enumerate(SeqIO.parse(src, "genbank"), 1):
        short = f"c{i}"
        rows.append(f"{short}\t{r.id}")
        r.id = r.name = short
        recs.append(r)
    SeqIO.write(recs, dest, "genbank")
    ids.write_text("\n".join(rows) + "\n")


def restore_json(path: Path, ids: Path) -> None:
    m = dict(line.split("\t") for line in ids.read_text().splitlines() if line)
    data = json.loads(path.read_text())
    for r in data.get("records", []):
        r["id"] = m.get(r["id"], r["id"])
    path.write_text(json.dumps(data))


if __name__ == "__main__":
    cmd, *args = sys.argv[1:]
    {"shorten": shorten, "restore-json": restore_json}[cmd](*map(Path, args))
