"""v3-active: the MIBiG 3.1-era entries (accession < BGC0002750) of objective 01's 5.0rc1 training archives, as a
MIBiG-format archive pair for ``sanntis train --external-gbk/-json``. GBK and JSON members are copied verbatim, so a
v3 BGC is the same locus in both objectives (incl. the 15 loci taken at their 5.0rc1 coordinates).

    python v3_archives.py <01 inputs/mibig_5.0rc1-train dir> <out dir>

writes <out>/mibig_v3-active_{gbk,json}.tar.gz, manifest.tsv (objective 01's rows), provenance.json
"""

import json
import sys
from pathlib import Path

import pandas as pd
from sanntis.train.collect import external_loci

from sv1.archives import read_tar, sha256, write_tar

V3_MAX = "BGC0002750"  # first accession after MIBiG 3.1
SRC = "mibig_5.0rc1-train"
NAME = "mibig_v3-active"


def main(src: Path, out: Path) -> None:
    man = pd.read_csv(src / "manifest.tsv", sep="\t", dtype={"locus_changed_since_4.0": bool})
    man = man[man["accession"] < V3_MAX].reset_index(drop=True)
    keep = set(man["accession"])
    arch = {}
    out.mkdir(parents=True, exist_ok=True)
    for kind in ("gbk", "json"):
        members = {n: b for n, b in read_tar(src / f"{SRC}_{kind}.tar.gz", f".{kind}").items()
                   if Path(n).stem in keep}
        if set(Path(n).stem for n in members) != keep:
            raise SystemExit(f"{kind}: {len(members)} members for {len(keep)} manifest rows")
        arch[kind] = out / f"{NAME}_{kind}.tar.gz"
        write_tar(arch[kind], members)
    man.to_csv(out / "manifest.tsv", sep="\t", index=False)

    n = len(external_loci(arch["gbk"], arch["json"], set()))  # what `sanntis train` reads, all-or-nothing
    assert n == len(man), (n, len(man))
    counts = {"packaged": len(man), "by_status": man["status"].value_counts().to_dict(),
              "by_gbk_source": man["gbk_source"].value_counts().to_dict(),
              "locus_changed_since_4.0": int(man["locus_changed_since_4.0"].sum()),
              "over_200kb": int((man["length"] > 200_000).sum())}
    prov = {"rule": f"accession < {V3_MAX}", "source": {p.name: sha256(p) for p in sorted(src.glob(f"{SRC}_*.tar.gz"))},
            "outputs": {p.name: sha256(p) for p in arch.values()}, "counts": counts}
    (out / "provenance.json").write_text(json.dumps(prov, indent=1) + "\n")
    print(json.dumps(counts))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
