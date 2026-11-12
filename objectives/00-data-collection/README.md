# 00 · Data collection

**Objective:** download every external input of objectives 01–03 from its public repository by accession,
verified against checksums.

```bash
bash objectives/00-data-collection/01_mibig.sh               # ~85 MB, minutes
bash objectives/00-data-collection/02_taxonomy.sh            # NCBI taxdump + GTDB r226 metadata, ~310 MB
bash objectives/00-data-collection/03_pool_genomes.sh        # negative pool v2 genomes
bash objectives/00-data-collection/04_testset_background.sh  # objective 02's test backgrounds
```

Trial run: `LIMIT=n DATA_ROOT=<empty dir> bash …` (the first n genomes per kingdom; MIBiG and the taxonomy always
in full). Parallel downloads: `JOBS=n`.

## Steps

| Script | Downloads | Verified by | Output (`$DATA_ROOT/raw/…`) |
|---|---|---|---|
| `01_mibig.sh` | MIBiG 4.0 JSON + GBK, 5.0rc1 JSON (dl.secondarymetabolites.org); 174 loci (173 distinct) of 5.0rc1 entries missing from the 4.0 GBK archive or whose locus changed since 4.0 (NCBI nuccore, E-utilities) | sha256 pinned in `mibig_archives.tsv` (the same pins as `sanntis.train.mibig`); nuccore: complete GenBank + count, `SHA256SUMS` written | `mibig/`, `mibig/nuccore/<acc>_<from>_<to>.gbk` |
| `02_taxonomy.sh` | NCBI taxdump `taxdmp_2026-07-01.zip` (NCBI taxdump archive); GTDB r226 `bac120` / `ar53` metadata | sha256 pinned in `taxonomy.tsv` | `ncbi/taxdump/`; the GTDB tables in `$SANNTIS_CACHE/gtdb/r226/`, where SanntiS reads them |
| `03_pool_genomes.sh` | the 897 genomes of negative pool v2 (`pool_v2_genomes.tsv`; NCBI Datasets): prokaryotes as genome FASTA, fungi with their annotated RefSeq GenBank | md5 manifest inside each NCBI package + count | `ncbi/assemblies/<acc>/<acc>.fna.gz` (+ `.gff.gz`, fungi: `.gbff.gz`, `assembly_data_report.jsonl`) |
| `04_testset_background.sh` | the 124 test-background genomes of objective 02 (`testset_background_genomes.tsv`; eukaryotes with their annotated GBFF) | as above | `ncbi/assemblies/<acc>/` |

A file that fails its checksum is deleted and the step stops; re-running downloads it again.

## Inputs (`accessions/`)

| File | Rows | Content |
|---|---|---|
| `mibig_archives.tsv` | 3 archives | SanntiS release pins |
| `mibig_nuccore_loci.tsv` | 174 | first locus of the 5.0rc1 active / pending entries whose GBK is not in the 4.0 archive (159) or whose first locus differs from 4.0 (15) |
| `taxonomy.tsv` | 3 files | NCBI taxdump and GTDB r226 metadata: URL, sha256, destination |
| `pool_candidates_pool_v1.tsv` | 217 (205 prokaryotes, 12 fungi) | candidates of the pool plan: the genomes of negative pool v1, `clean_cds` = their BGC-free CDS in pool v1, priority 0 |
| `pool_candidates_bacdive.tsv` | 2,521 | candidates: BacDive type strains (BacDive 2026-03-23) with an NCBI assembly, `clean_cds` = 0.9 x their CDS count, priority 0 |
| `pool_candidates_refseq_fungi.tsv` | 674 | candidates: RefSeq fungal assemblies (`assembly_summary.txt` of 2026-10-29), `clean_cds` = genome size / 2,900 x 0.95, priority 1 |
| `pool_breadth.tsv` | 205 | the pool v1 prokaryotes; those not planned join the pool with at most 500 CDS each |
| `pool_exclude.tsv` | 254 | genomes never in the pool: 129 species held out for genome-level evaluation, the 124 test-background genomes (assemblies; eukaryotes also their species), and one GTDB r226 representative suppressed by NCBI since the release (GCA_946479335.1) |
| `pool_v2_genomes.tsv` | 897 (735 prokaryotes, 162 fungi) | the plan of negative pool v2 (`sanntis plan-negatives`, objective 01 `02_plan.sh`): genome, accession, lineage, `quota_cds`, plan columns |
| `testset_background_genomes.tsv` | 124 (102 GTDB r226 prokaryotes, 22 RefSeq fungi / protozoa) | objective 02 `workers/sample_background.py`: no assembly and no species of the pool |

Every candidate table carries the NCBI taxid; lineages come from the pinned taxdump inside `sanntis
plan-negatives`. GTDB r226 species representatives (CheckM2 completeness >= 90 %, contamination <= 5 %, priority 1)
are added by `plan-negatives` from the GTDB metadata.

## Notes

- **MIBiG nuccore loci.** 159 loci are missing from the 4.0 archive; the 15 loci changed since 4.0 (5 only in
  notation, `0..0` vs `1..<end>`) are taken at their 5.0rc1 coordinates. BGC0003182 and BGC0003249 have no row:
  neither is in the 4.0 archive and both drop out of training (WGS contig with no CDS; 8.4 Mb whole genome).
- **Pool plan and test backgrounds.** The plan is computed first, with the held-out species excluded; the test
  backgrounds are then sampled with no assembly (and no GTDB species, for prokaryotes; no NCBI species, for
  eukaryotes) of the planned pool, and added to `pool_exclude.tsv`. The plan's order within a lineage group depends
  only on the seed and the genome, so the plan is the same with the final exclusion table (objective 01 checks it).
- **Snapshots.** The BacDive, RefSeq and GTDB candidate lists and the background draw are inputs of record: the
  RefSeq `assembly_summary.txt` changes daily and is not archived by NCBI, so the stored tables are what a rebuild
  uses.
