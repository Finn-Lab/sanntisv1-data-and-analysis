# 00 · Data collection

**Objective:** download every external input of objectives 01–03 from its public repository by accession,
verified against checksums.

```bash
bash objectives/00-data-collection/01_mibig.sh               # ~85 MB, minutes
bash objectives/00-data-collection/02_pool_genomes.sh        # optional (pool v1 itself is on Zenodo)
bash objectives/00-data-collection/03_testset_backbones.sh   # optional (the built test sets are on Zenodo)
bash objectives/00-data-collection/05_testset_background.sh   # objective 02's test backgrounds
```

Trial run on the first rows of every table: `LIMIT=3 DATA_ROOT=<scratch dir> bash …`. Parallel downloads: `JOBS=n`.

## Steps

| Script | Downloads | Verified by | Output (`$DATA_ROOT/raw/…`) |
|---|---|---|---|
| `01_mibig.sh` | MIBiG 4.0 JSON + GBK, 5.0rc1 JSON (dl.secondarymetabolites.org); 174 loci (173 distinct) of 5.0rc1 entries missing from the 4.0 GBK archive or whose locus changed since 4.0 (NCBI nuccore, E-utilities) | sha256 pinned in `mibig_archives.tsv` (the same pins as `sanntis.train.mibig`); nuccore: complete GenBank + count, `SHA256SUMS` written | `mibig/`, `mibig/nuccore/<acc>_<from>_<to>.gbk` |
| `02_pool_genomes.sh` | 228 genomes of negative pool v1 (NCBI Datasets) | md5 manifest inside each NCBI package + count | `ncbi/assemblies/<acc>/<acc>.fna.gz` (+ `.gff.gz`, `assembly_data_report.jsonl`) |
| `03_testset_backbones.sh` | 138 backbone genomes of the nucleotide test sets (NCBI Datasets) | as above | `ncbi/assemblies/<acc>/` |
| `05_testset_background.sh` | 124 test-background genomes of objective 02, none in pool v1 (NCBI Datasets; eukaryotes with their annotated GBFF) | as above | `ncbi/assemblies/<acc>/` (+ `.gbff.gz`) |

A file that fails its checksum is deleted and the step stops; re-running downloads it again.

## Inputs (`accessions/`)

| File | Rows | Source of the list |
|---|---|---|
| `mibig_archives.tsv` | 3 archives | SanntiS release pins |
| `mibig_nuccore_loci.tsv` | 174 | first locus of the 5.0rc1 active / pending entries whose GBK is not in the 4.0 archive (159) or whose first locus differs from 4.0 (15) |
| `pool_v1_genomes.tsv` | 228 (205 GTDB r226 proka, 23 RefSeq euk) | negative pool v1 manifest |
| `testsets_nucleotide_backbones.tsv` | 138 assemblies + source contigs | v5 isolate (55 + 15) and v5shuf2 (102 + 22) builds (all pool v1 genomes) |
| `testset_background_genomes.tsv` | 124 (102 GTDB r226 proka, 22 RefSeq fungi / protozoa) | objective 02 `workers/sample_background.py` (2026-10-29): pool v1 removed (assemblies; eukaryotes also its species), then the build-negatives sampling |

## Notes

- **MIBiG nuccore loci.** The 159 loci missing from the 4.0 archive equal what `sanntis train` fetched for the release
  run (same sequences; SanntiS then dropped CDS cut by the locus slice and pseudo genes, so its cached copies differ
  only in those features). The 15 loci changed since 4.0 (5 only in notation, `0..0` vs `1..<end>`) were added for
  objective 01, which takes them at their 5.0rc1 coordinates; the release run used the 4.0 GBKs for them.
  BGC0003182 and BGC0003249 have no row: neither is in the 4.0 archive and both drop out of training (WGS contig with
  no CDS; 8.4 Mb whole genome).
- **Test-background genomes** share no assembly (and, for eukaryotes, no species) with pool v1; GTDB representatives
  are one per species, so neither do the prokaryotes. 11 of the 102 GTDB genera also occur in pool v1. The RefSeq
  draw used the 2026-10-29 `assembly_summary.txt`; the stored list is the reproducible input.
