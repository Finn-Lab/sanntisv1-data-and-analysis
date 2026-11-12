# 01 · Production models

**Objective:** build negative pool v2 and train the SanntiS 1.2.0 production bundle with `sanntis train` on the
MIBiG 5.0rc1 active + pending entries (GBKs from the MIBiG 4.0 archive, else NCBI nuccore) and pool v2, with the
default protocol (CV + tuning on the training pool itself).

```bash
bash objectives/01-production-models/01_inputs.sh   # MIBiG-format archives, ~2 min
bash objectives/01-production-models/02_plan.sh     # InterProScan 6 on the positives + the genome plan of pool v2
bash objectives/01-production-models/03_pool.sh     # sanntis build-negatives → pool v2 (+ the Hub deposit copy)
bash objectives/01-production-models/04_train.sh    # sanntis train (NEGATIVES=<pool dir | Hub revision>)
```

Status: TODO (filled from the run).

## Inputs

MIBiG 5.0rc1 JSON + 4.0 JSON / GBK archives and the nuccore loci, the pool genomes, the candidate and exclusion
tables, the NCBI taxdump and the GTDB r226 metadata (objective 00); InterProScan 6.0.3.0 / InterPro 107.0 data
(downloaded by SanntiS on first use). NPC labels are computed by `sanntis train` from the compounds of the JSON
entries (NPClassifier API, cached in `SANNTIS_CACHE`).

## Steps

| Script | Does | Output (`$DATA_ROOT/01/…`) |
|---|---|---|
| `01_inputs.sh` | worker `mibig_train_archives.py`: 5.0rc1 active + pending as a MIBiG-format archive pair, checked with SanntiS's reader for external data and against `expected/`; packs the Zenodo file | `inputs/mibig_5.0rc1-train/`, `deposit/sanntisv1_mibig5.0rc1_train.tar.gz` |
| `02_plan.sh` | `sanntis train --until annotate` (the deduplicated production positives and their InterProScan 6 tokens), then `sanntis plan-negatives` (below); the plan must equal objective 00's `pool_v2_genomes.tsv` | `plan/{run/, genomes.tsv, genomes.stats.json}` |
| `03_pool.sh` | `sanntis build-negatives --genomes <plan> --exclude pool_exclude.tsv --pool-version v2` on the downloaded genomes; the Hub deposit copy | `pool/{prepared/, antismash/, ips6/, v2/}`, `hub/sanntis-training-data/negatives/v2/` |
| `04_train.sh` | `sanntis train --mibig none --external-gbk … --external-json … --negatives pool/v2`; checks the collected pool (2,777 loci) | `run/` (`bundle/`, `report/report.md`, `tune/tune.json`) |

## Training data

### Positives (`01_inputs.sh`)

MIBiG 5.0rc1 ships JSON only, so `sanntis train --mibig 5.0rc1` has no GBK archive to read. The 5.0rc1 entries go
in as external data in the MIBiG format (`mibig_5.0rc1-train_{gbk,json}.tar.gz`, one `<accession>.gbk` /
`<accession>.json` per BGC):

| | Entries | |
|---|---|---|
| 5.0rc1 active + pending | 2,797 (2,598 + 199) | pending = active in 4.0, under re-review in 5.0 |
| GBK from the 4.0 archive, verbatim | 2,621 | first locus unchanged since 4.0 |
| GBK from nuccore, first 5.0rc1 locus | 174 | 159 not in the 4.0 archive + 15 whose locus changed since 4.0; `/pseudo` CDS, CDS cut by the slice and repeated CDS locations dropped |
| excluded (`excluded.tsv`) | 2 | BGC0003182 (WGS contig without CDS), BGC0003249 (8.4 Mb genome) |
| **packaged** | **2,795** | `sanntis train` then excludes the 18 loci > 200 kb → **2,777 loci** |

The archives are written deterministically. `expected/mibig_train_manifest.tsv` lists every BGC (status, GBK source,
locus, length, CDS, md5 of the sequence and of the proteins as SanntiS reads them) and is what a rebuild must
reproduce; `env/zenodo.tsv` (`mibig5_train`) pins the deposit file, reproduced by a build with the same SanntiS and
Biopython (`provenance.json` records both).

### Negative pool v2 (`02_plan.sh`, `03_pool.sh`)

The pool follows the lineages of the positives: per lineage group, negative CDS in proportion to the positives' CDS.

- **Plan** (`sanntis plan-negatives`). Positives: the production positives of the plan run (deduplicated; fungi
  only among eukaryotes), lineage from their MIBiG NCBI taxid in the NCBI taxdump of 2026-07-01. Budgets:
  469,554 prokaryotic CDS and 186,849 fungal CDS (41.12 x the fungal positives' CDS). One genome with 250 CDS for
  each phylum with candidates but no positives; the rest per lineage group (the positives' family when it has
  candidates, else their order, else their phylum), at most 1,500 CDS and one genome per species, quotas >= 100 CDS,
  seed 0. Candidates: pool v1 genomes and BacDive type strains (priority 0), GTDB r226 species representatives
  (CheckM2 >= 90 % / <= 5 %) and RefSeq fungi (priority 1); `pool_exclude.tsv` removed first. Every pool v1
  prokaryote not planned joins with at most 500 CDS.
- **Build** (`sanntis build-negatives`). Prokaryotic genes called with pyrodigal (single mode; meta below 100 kb),
  fungal CDS from the RefSeq GenBank (one CDS per gene); antiSMASH 8.0.1 `--minimal --genefinding-tool none`
  (`--taxon bacteria` / `fungi`); CDS in an antiSMASH region removed, the rest forms BGC-free runs; the longest run of
  >= 50 CDS of 55 prokaryotic and 15 fungal genomes is held out as an isolate test backbone (objective 02's isolate
  set); per genome, `quota_cds` BGC-free CDS drawn as 3 portions of the concatenated runs; InterProScan 6.0.3.0 /
  InterPro 107.0 (14 member databases) on the portion and backbone proteins.

| | Genomes | Quota CDS |
|---|---|---|
| prokaryotes | 735 (367 GTDB r226, 163 BacDive, 205 pool v1) | 517,515 |
| fungi | 162 (157 RefSeq, 5 pool v1) | 184,731 |

Counts after the build (CDS, portions, backbones, InterProScan coverage): TODO (from the run).

### Recipe

SanntiS 1.2.0's recipe (`recipe.yaml`): fungi only among eukaryotic positives and pool genomes; 30 % of the
positives placed in runs of 2–4 BGCs separated by 2–40 background CDS; decoding by the mean ORF score, CDS without
InterPro hits scored as the lower of their annotated neighbours, hysteresis gap 0.05, no merge gap; operating point
t* = 0.20 (low threshold 0.15).

## Expected

`expected/mibig_train_{manifest,excluded}.tsv`: the training archives (above). The bundle report tables of the run
(`tune.tsv`, `validation_oof.tsv`, `validation.tsv` with 95 % CIs, strata, class curves, sweep, NPC validation,
collected-pool counts, bundle manifest / config): TODO (frozen from the run).
