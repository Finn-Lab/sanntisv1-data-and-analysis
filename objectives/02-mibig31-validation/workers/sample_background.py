"""Test-background genomes of the shuffled test set, disjoint from the training negative pool v2: run once by the
maintainers after the pool plan; the selection is frozen as
``objectives/00-data-collection/accessions/testset_background_genomes.tsv`` (the reproducible input) and its
assemblies are added to ``pool_exclude.tsv``.

Same sampling code as ``sanntis build-negatives`` (``sanntis.negpool.sample``), with the pool removed first:
- prokaryotes: GTDB r226 species representatives of no GTDB species of a pool assembly (any GCA / GCF version),
  ``balanced_sample`` over phyla, ``n_proka`` genomes, seed 43;
- eukaryotes: RefSeq fungi + protozoa not in the pool and of none of its NCBI species, ``REFSEQ_GROUPS`` genomes
  per phylum (seed 1), then ``n_euk`` of them drawn per group in proportion.

    python sample_background.py <pool_v2_genomes.tsv> <out tsv> [n_proka=102] [n_euk=22]
"""

import re
import sys
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
from sanntis.negpool import sample as S
from sanntis.train.ncbi import LineageCache
from sanntis.utils import cache_dir

COLUMNS = ["accession", "genome", "kingdom", "group", "source", "taxon", "lineage"]


def core(acc: str) -> str:
    """GCA_000220375.1 / RS_GCF_000220375.2 → 000220375 (the assembly, whatever the version or database)."""
    m = re.search(r"GC[AF]_(\d+)", acc)
    return m.group(1) if m else acc


def proka(pool: pd.DataFrame, n: int) -> pd.DataFrame:
    df = pd.concat([pd.read_csv(t, sep="\t", dtype=str, low_memory=False) for t in S.gtdb_tables("226")],
                   ignore_index=True)
    taken = {core(a) for a in pool["accession"]}
    # GTDB species of the pool's assemblies (representatives or not): no background shares one
    species = set(df.loc[df["accession"].map(core).isin(taken), "gtdb_taxonomy"].str.extract(r"s__([^;]+)")[0]
                  .dropna())
    df = df[df["gtdb_representative"].str.lower().eq("t")].copy()
    df = df[~df["accession"].map(core).isin(taken)
            & ~df["gtdb_taxonomy"].str.extract(r"s__([^;]+)")[0].isin(species)]
    df = df.join(S.split_gtdb_taxonomy(df["gtdb_taxonomy"]))
    df = df[df["phylum"].notna()]
    s = S.balanced_sample(df, "phylum", n, seed=43)
    acc = s["accession"].str.replace(r"^(RS_|GB_)", "", regex=True)
    return pd.DataFrame({"accession": acc, "genome": acc, "kingdom": "proka", "group": s["domain"].to_numpy(),
                         "source": "gtdb", "taxon": "bacteria", "lineage": s["gtdb_taxonomy"].to_numpy()})


def euk(pool: pd.DataFrame, n: int) -> pd.DataFrame:
    groups = ["fungi", "protozoa"]
    root = cache_dir("refseq") / datetime.now(UTC).strftime("%Y-%m-%d")
    summaries = {g: S.fetch(S.REFSEQ_SUMMARY.format(group=S.REFSEQ_GROUPS[g][1]), root / f"assembly_summary_{g}.txt")
                 for g in groups}
    taken = {core(a) for a in pool["accession"]}
    species = set(pool.loc[pool["kingdom"] == "euk", "species"]) - {""}
    trimmed = {}
    for g, p in summaries.items():  # the pool out before sampling, so the per-phylum draw fills from the rest
        d = S.load_assembly_summary(p)
        d = d[~d["assembly_accession"].map(core).isin(taken) & ~d["organism_name"].fillna("").map(
            lambda o: any(o == s or o.startswith(s + " ") for s in species))]
        trimmed[g] = root / f"assembly_summary_{g}.no_pool.txt"
        head = p.read_text().splitlines(keepends=True)[:2]  # the reader expects the two header lines
        with open(trimmed[g], "w") as fh:
            fh.writelines(head)
            d.to_csv(fh, sep="\t", index=False, header=False)
    taxids = {int(t) for p in trimmed.values() for t in S.load_assembly_summary(p)["taxid"].dropna()
              if str(t).isdigit()}
    s = S.sample_refseq(trimmed, LineageCache().resolve(taxids), seed=1)
    if len(s) > n:  # proportional per group, seeded
        quota = (s["group"].value_counts() / len(s) * n).round().astype(int)
        quota.iloc[0] += n - quota.sum()
        s = pd.concat([s[s["group"] == g].sample(n=int(q), random_state=1) for g, q in quota.items()])
    return s.assign(genome=s["accession"])[COLUMNS]


def main(pool_tsv: Path, out: Path, n_proka: int, n_euk: int) -> None:
    pool = pd.read_csv(pool_tsv, sep="\t", dtype=str, keep_default_na=False)
    df = pd.concat([proka(pool, n_proka), euk(pool, n_euk)], ignore_index=True)
    df = df.sort_values(["kingdom", "accession"], ascending=[False, True])
    assert not set(df["accession"].map(core)) & {core(a) for a in pool["accession"]}
    df[COLUMNS].to_csv(out, sep="\t", index=False)
    print(df.groupby(["kingdom", "group"]).size().to_string())


if __name__ == "__main__":
    a = sys.argv[1:]
    main(Path(a[0]), Path(a[1]), int(a[2]) if len(a) > 2 else 102, int(a[3]) if len(a) > 3 else 22)
