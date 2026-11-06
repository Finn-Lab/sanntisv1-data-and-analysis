"""MIBiG 5.0rc1 active + pending entries as a MIBiG-format archive pair for ``sanntis train --external-gbk/-json``.

MIBiG 5.0rc1 ships JSON only, and ``sanntis train --mibig <release>`` needs a GBK archive of the release, so the
5.0rc1 training data goes in as external data. The GBK of an entry is

- the 4.0 archive GBK, verbatim, when its first locus is unchanged since 4.0;
- otherwise the first 5.0rc1 locus from NCBI nuccore (objective 00), cleaned as the release run did (``/pseudo`` CDS,
  CDS cut by the locus slice and repeated CDS locations dropped; the record renamed to the accession).

JSON entries are the 5.0rc1 files, verbatim. Entries without a usable GBK are excluded (``excluded.tsv``); loci
> 200 kb stay in, ``sanntis train`` excludes them itself. The archives are checked with SanntiS's own reader for
external data (``sanntis.train.collect.external_loci``) and written deterministically (sorted members, mtime 0).

    python mibig_train_archives.py <00 raw/mibig dir> <out dir>

writes <out>/mibig_5.0rc1-train_{gbk,json}.tar.gz, manifest.tsv, excluded.tsv, provenance.json
"""

import hashlib
import io
import json
import sys
from collections import Counter
from importlib.metadata import version
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from sanntis.io.genbank import iter_cds
from sanntis.train import mibig
from sanntis.train.collect import external_loci

from sv1.archives import read_tar, sha256, write_tar

STATUSES = ("active", "pending")
NAME = "mibig_5.0rc1-train"


def first_locus(entry: dict) -> tuple[str, int, int] | None:
    loci = entry.get("loci") or []
    if not loci:
        return None
    x = loci[0]
    return x["accession"], int(x["location"]["from"]), int(x["location"]["to"])


def clean_nuccore(rec: SeqRecord, accession: str) -> SeqRecord:
    """As the MIBiG GBKs (the release run's ``clean_nuccore_record``): drop /pseudo CDS, CDS cut by the locus slice
    (bases cannot encode the /translation) and repeated CDS locations (alternative transcripts; first kept)."""
    seen: set[str] = set()
    kept = []
    for f in rec.features:
        if f.type == "CDS":
            if "pseudo" in f.qualifiers:
                continue
            tr = (f.qualifiers.get("translation") or [""])[0]
            if tr and len(f.location) < 3 * len(tr):
                continue
            if str(f.location) in seen:
                continue
            seen.add(str(f.location))
        kept.append(f)
    rec.features = kept
    rec.id = rec.name = accession
    rec.annotations["accessions"] = [accession]
    rec.annotations.setdefault("molecule_type", "DNA")
    return rec


def main(raw: Path, out: Path) -> None:
    src = {"json_5.0rc1": raw / "mibig_json_5.0rc1.tar.bz2", "json_4.0": raw / "mibig_json_4.0.tar.gz",
           "gbk_4.0": raw / "mibig_gbk_4.0.tar.gz", "nuccore_sums": raw / "nuccore" / "SHA256SUMS"}
    v5 = mibig.read_entries(src["json_5.0rc1"])
    v4 = mibig.read_entries(src["json_4.0"])
    sel = {a: e for a, e in v5.items() if e.get("status") in STATUSES}
    gbk4 = mibig.read_gbk_archive(src["gbk_4.0"], set(sel))
    raw_json = {json.loads(b)["accession"]: b for b in read_tar(src["json_5.0rc1"], ".json").values()}  # as published

    gbks, rows, excl = {}, [], []
    for acc in sorted(sel):
        e, loc = sel[acc], first_locus(sel[acc])
        changed = acc in v4 and first_locus(v4[acc]) != loc
        row = {"accession": acc, "status": e["status"], "locus": "{}:{}-{}".format(*loc) if loc else "",
               "locus_changed_since_4.0": changed}
        try:
            if errs := mibig.validate_entry(e, mibig=True):
                raise ValueError("json: " + "; ".join(errs))
            if acc in gbk4 and not changed:
                row["gbk_source"], data = "mibig_gbk_4.0", gbk4[acc]
                rec = mibig.archive_record(data)
            else:
                f = raw / "nuccore" / ("{}_{}_{}.gbk".format(*loc) if loc else "-")
                if not f.is_file():
                    raise ValueError("no GBK: not in the 4.0 archive and no objective-00 nuccore locus")
                rec = clean_nuccore(SeqIO.read(f, "genbank"), acc)
                buf = io.StringIO()
                SeqIO.write(rec, buf, "genbank")
                row["gbk_source"], data = "nuccore", buf.getvalue().encode()
                rec = mibig.archive_record(data)
            cds = list(iter_cds(rec, missing_translation="error"))
            if not cds:
                raise ValueError(f"no CDS ({row['gbk_source']})")
        except Exception as x:
            excl.append({"accession": acc, "status": e["status"], "reason": str(x)[:200]})
            continue
        gbks[f"{acc}.gbk"] = data
        rows.append({**row, "length": len(rec.seq), "n_cds": len(cds),
                     "seq_md5": hashlib.md5(str(rec.seq).upper().encode()).hexdigest(),
                     "proteins_md5": hashlib.md5("\n".join(c.md5 for c in cds).encode()).hexdigest()})

    out.mkdir(parents=True, exist_ok=True)
    arch = {"gbk": out / f"{NAME}_gbk.tar.gz", "json": out / f"{NAME}_json.tar.gz"}
    write_tar(arch["gbk"], gbks)
    write_tar(arch["json"], {f"{r['accession']}.json": raw_json[r["accession"]] for r in rows})
    man = pd.DataFrame(rows)[["accession", "status", "gbk_source", "locus", "locus_changed_since_4.0", "length",
                              "n_cds", "seq_md5", "proteins_md5"]]
    man.to_csv(out / "manifest.tsv", sep="\t", index=False)
    pd.DataFrame(excl, columns=["accession", "status", "reason"]).to_csv(out / "excluded.tsv", sep="\t", index=False)

    # the gate: what `sanntis train --external-gbk/--external-json` will read, all-or-nothing
    n = len(external_loci(arch["gbk"], arch["json"], set()))
    assert n == len(rows), (n, len(rows))

    counts = {"selected": dict(Counter(e["status"] for e in sel.values())), "packaged": len(rows),
              "excluded": len(excl), "by_status": man["status"].value_counts().to_dict(),
              "by_gbk_source": man["gbk_source"].value_counts().to_dict(),
              "locus_changed_since_4.0": int(man["locus_changed_since_4.0"].sum()),
              "over_200kb": int((man["length"] > 200_000).sum())}
    prov = {"statuses": list(STATUSES), "inputs": {k: {"name": p.name, "sha256": sha256(p)} for k, p in src.items()},
            "outputs": {p.name: sha256(p) for p in arch.values()}, "counts": counts,
            "versions": {"sanntis": version("sanntis"), "biopython": version("biopython")}}
    (out / "provenance.json").write_text(json.dumps(prov, indent=1) + "\n")
    print(json.dumps(counts))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
