"""The test-background pool: objective 00's test-background genomes through the stages of ``sanntis build-negatives``
(its own units and code) up to the antiSMASH stripping; no InterProScan (the test contigs are nucleotide GenBank,
SanntiS annotates their proteins at scoring time).

- clean (eukaryotes, here): the RefSeq GBFF without organelle records (mitochondrion, plastid, ...), keeping only
  source + CDS features; a CDS repeating an earlier one's location or name (locus_tag, else gene, else protein_id;
  alternative transcripts) or with overlapping exons is dropped. antiSMASH 8.0.1 rejects such records; counts per
  genome in ``cleaned/cleaning.tsv``. A genome left with no CDS (an organelle-only assembly) is not a background
  (``cleaned/dropped.tsv``)
- prepare (label ``cpu_small``): prokaryotes: pyrodigal on the NCBI genome (single mode; meta below 100 kb);
  eukaryotes: the cleaned GBFF CDS → ``prepared/<genome>/{input.gbk, cds.parquet, proteins.faa}``
- antiSMASH 8.0.1 (label ``antismash``): ``--minimal --genefinding-tool none``, taxon bacteria / fungi
- strip: CDS in antiSMASH regions → native mask; the rest → BGC-free runs; the isolate test backbones (longest run
  >= 50 CDS of 55 prokaryotic + 15 eukaryotic genomes, the pool's §2.5 rule)

    python background_pool.py <genome list> <assemblies dir> <build dir> [--runner nextflow|inline] [--workers N]
        [--container-engine E] [--nf-config F ...] [--antismash-databases D]

writes <build dir>/background/{genomes,cds,clean_runs,native_mask,test_backbones}.parquet + manifest.json
"""

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from sanntis.config import NextflowOptions
from sanntis.io.genbank import read_records
from sanntis.negpool import build as B
from sanntis.negpool import sample as S
from sanntis.runner.tasks import RunnerOptions, Unit, dispatch, write_unit_script

from sv1.archives import sha256


def overlapping_exons(loc) -> bool:
    parts = sorted(getattr(loc, "parts", [loc]), key=lambda p: int(p.start))
    return any(int(b.start) < int(a.end) for a, b in zip(parts, parts[1:], strict=False))


def clean_gbff(src: Path, dest: Path) -> dict:
    """RefSeq GBFF → records antiSMASH 8.0.1 accepts (see module doc); returns what was dropped."""
    n = {"organelle_records": 0, "repeated_location_cds": 0, "repeated_name_cds": 0, "overlapping_exon_cds": 0,
         "cds": 0}
    recs = []
    for rec in read_records(src):
        src_feat = next((f for f in rec.features if f.type == "source"), None)
        if src_feat is not None and "organelle" in src_feat.qualifiers:
            n["organelle_records"] += 1
            continue
        seen, names, keep = set(), set(), [src_feat] if src_feat is not None else []
        for f in rec.features:
            if f.type != "CDS":
                continue
            name = next((f.qualifiers[k][0] for k in ("locus_tag", "gene", "protein_id") if f.qualifiers.get(k)), "")
            if overlapping_exons(f.location):
                n["overlapping_exon_cds"] += 1
            elif str(f.location) in seen:
                n["repeated_location_cds"] += 1
            elif name and name in names:  # another isoform of a gene already kept
                n["repeated_name_cds"] += 1
            else:
                seen.add(str(f.location))
                names.add(name)
                keep.append(f)
        rec.features = keep
        n["cds"] += len(keep) - (src_feat is not None)
        recs.append(rec)
    tmp = dest.with_suffix(".part")
    with gzip.open(tmp, "wt") as fh:
        SeqIO.write(recs, fh, "genbank")
    tmp.replace(dest)
    return n


def genome_table(listing: Path, assemblies: Path, out: Path) -> pd.DataFrame:
    g = pd.read_csv(listing, sep="\t", dtype=str, keep_default_na=False)
    ext = g["kingdom"].map({"proka": "fna.gz", "euk": "gbff.gz"})
    g["path"] = [str(assemblies / a / f"{a}.{e}") for a, e in zip(g["accession"], ext, strict=True)]
    missing = [p for p in g["path"] if not Path(p).is_file()]
    if missing:
        raise SystemExit(f"{len(missing)} genome files missing (run objectives/00-data-collection/"
                         f"05_testset_background.sh), e.g. {missing[0]}")
    cdir = out / "cleaned"
    cdir.mkdir(exist_ok=True)
    rows = []
    for i in g.index[g["kingdom"] == "euk"]:
        dest = cdir / f"{g.at[i, 'genome']}.gbff.gz"
        if not dest.exists():
            rows.append({"genome": g.at[i, "genome"], **clean_gbff(Path(g.at[i, "path"]), dest)})
        g.at[i, "path"] = str(dest)
    if rows:
        log = cdir / "cleaning.tsv"
        old = pd.read_csv(log, sep="\t") if log.exists() else None
        pd.concat([old, pd.DataFrame(rows)]).drop_duplicates("genome", keep="last").sort_values("genome").to_csv(
            log, sep="\t", index=False)
    g["url"] = ""
    # assemblies with no nuclear CDS (e.g. RefSeq organelle-only genomes) cannot be a background: listed, left out
    clean = pd.read_csv(cdir / "cleaning.tsv", sep="\t") if (cdir / "cleaning.tsv").exists() else pd.DataFrame()
    empty = set(clean.loc[clean["cds"] == 0, "genome"]) if len(clean) else set()
    if empty:
        g[g["genome"].isin(empty)].to_csv(cdir / "dropped.tsv", sep="\t", index=False)
        print(f"left out, no nuclear CDS: {sorted(empty)}")
        g = g[~g["genome"].isin(empty)].reset_index(drop=True)
    return g.reindex(columns=S.GENOME_COLUMNS, fill_value="")


def main(a: argparse.Namespace) -> None:
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    genomes = genome_table(Path(a.genomes), Path(a.assemblies), out)
    genomes.to_csv(out / "genomes.tsv", sep="\t", index=False)  # read by build-negatives' units
    nf = NextflowOptions(container_engine=a.container_engine, nf_config=a.nf_config or [])
    ropts = RunnerOptions(runner=a.runner, workers=a.workers, nextflow=nf)
    dispatch(out, "negpool_prepare", [Unit(g, "cpu_small") for g in genomes["genome"]], ropts)

    if not a.container_engine:
        B.check_host_antismash()
    db = Path(a.antismash_databases) if a.antismash_databases else None
    for g, taxon in zip(genomes["genome"], genomes["taxon"], strict=True):
        write_unit_script(out, "negpool_antismash", g, B.antismash_script(out, g, taxon, db))
    ropts.binds = [str(db)] if db else []
    ropts.container_prefix = B.container_prefix(a.container_engine, None, [out, *([db] if db else [])])
    dispatch(out, "negpool_antismash", [Unit(g, "antismash") for g in genomes["genome"]], ropts)

    tables = B.pool_tables(out, genomes)
    pool = out / "background"
    pool.mkdir(exist_ok=True)
    frames = {"genomes.parquet": genomes[["genome", "kingdom", "lineage", "source", "accession", "group"]],
              "cds.parquet": tables.pop("cds"), **tables}
    for name, df in frames.items():
        df.reset_index(drop=True).to_parquet(pool / name, index=False)
    kingdom = dict(zip(genomes["genome"], genomes["kingdom"], strict=True))
    counts = {k: {"genomes": int((genomes["kingdom"] == k).sum()),
                  "cds": int((frames["cds.parquet"]["genome"].map(kingdom) == k).sum()),
                  "clean_runs": int((frames["clean_runs.parquet"]["genome"].map(kingdom) == k).sum()),
                  "native_cds": int((frames["native_mask.parquet"]["genome"].map(kingdom) == k).sum()),
                  "test_backbones": int((frames["test_backbones.parquet"]["kingdom"] == k).sum())}
              for k in ("proka", "euk")}
    ver = sorted({(out / "antismash" / g / "version.txt").read_text().strip() for g in genomes["genome"]})
    man = {"antismash": {"version": ver, "options": "--minimal --genefinding-tool none --taxon bacteria|fungi"},
           "gene_calling": {"proka": "pyrodigal (single mode; meta below 100 kb)", "euk": "RefSeq GBFF CDS"},
           "euk_cleaning": pd.read_csv(out / "cleaned" / "cleaning.tsv", sep="\t").drop(columns="genome").sum()
           .astype(int).to_dict(),
           "counts": counts, "files": {n: sha256(pool / n) for n in frames}}
    (pool / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    print(json.dumps(counts))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("genomes")
    p.add_argument("assemblies")
    p.add_argument("out")
    p.add_argument("--runner", default="nextflow", choices=["nextflow", "inline"])
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--container-engine")
    p.add_argument("--nf-config", action="append")
    p.add_argument("--antismash-databases")
    main(p.parse_args())
