# 02 · MIBiG 3.1 validation

**Objective:** train on v3-active (MIBiG 3.1-era entries, accession < BGC0002750, active / pending in 5.0rc1: 2,419
BGCs) with negative pool v2 and score the model on Test A (5.0rc1 active, new since v3: 178) and Test B (active +
pending: 375), each inserted into isolate backbones held out of the pool and into gene-shuffled contigs of genomes
**that are not in the pool**; strata by class, lineage and novelty.

```bash
bash objectives/02-mibig31-validation/01_v3active.sh     # v3-active MIBiG-format archives (from objective 01's)
bash objectives/02-mibig31-validation/02_train.sh        # sanntis train → the v3-active model
bash objectives/02-mibig31-validation/03_background.sh   # test-background pool (gene calls, antiSMASH 8.0.1, strip)
bash objectives/02-mibig31-validation/04_testsets.sh     # isolate + shuffled nucleotide test sets, Zenodo deposit
bash objectives/02-mibig31-validation/05_score.sh        # InterProScan 6 tokens of the test proteins, Test A / B scores
bash objectives/02-mibig31-validation/06_report.sh       # report.md + figures, the v3active_eval_v3 deposit
```

Status: TODO (filled from the run).

## Steps

| Script | Does | Output (`$DATA_ROOT/02/…`) |
|---|---|---|
| `01_v3active.sh` | worker `v3_archives.py`: the < BGC0002750 members of objective 01's archives, verbatim; checked with SanntiS's external-data reader and against `expected/v3active_manifest.tsv` | `inputs/mibig_v3-active/` |
| `02_train.sh` | `sanntis train --mibig none --external-gbk/-json` on v3-active + pool v2 (objective 01); checks the collected pool (2,402 loci) | `run/` (`bundle/`, `report/`, `tune/`) |
| `03_background.sh` | worker `background_pool.py`: objective 00's 124 test-background genomes through `sanntis build-negatives`' stages up to the stripping; no InterProScan | `background/background/*.parquet` |
| `04_testsets.sh` | worker `build_testsets.py`: the Test B BGCs inserted into pool v2's isolate test backbones and into gene-shuffled background contigs; the deposit tarball (`env/zenodo.tsv`: `nucleotide_testsets_v3`) | `testsets/{testsets/, test_bgcs.parquet, contigs.tsv, manifest.json}`, `deposit/` |
| `05_score.sh` | step `02_tokens`, worker `tokens.py`: InterProScan 6 through SanntiS's annotation (md5 cache) → tokens + Pfam sets per test protein; step `02_score`, worker `score_testsets.py`: Test A / B on both sets at the bundle's t* with `sanntis.train.metrics` (CIs, strata, radar, PR curves, NPC), novelty vs the training positives, the test BGCs' NPC labels from objective 01's `collect` | `eval/{tokens.parquet, test_bgcs.parquet, labels/, bundle → run/bundle, score/{A,B}-{isolate,shuffled}.pkl, score/test_*.tsv}` |
| `06_report.sh` | step `02_report`, worker `report_v3active.py`: training and OOF-validation tables, Test A / B metrics with CIs, threshold sweeps, strata radars, novelty, NPC; step `02_deposit`: the `v3active_eval_v3` deposit (bundle, test-BGC table, labels, tokens, scores) | `report/{report.md, summary.tsv, figures/}`, `deposit/sanntisv1_v3active_eval_v3.tar.gz` |

## Training data (`01_v3active.sh`)

| | BGCs |
|---|---|
| objective 01 archives (5.0rc1 active + pending) with accession < BGC0002750 | 2,419 (2,418 active, 1 pending) |
| of which loci from nuccore | 18 (incl. the 15 at their changed 5.0rc1 coordinates) |
| `sanntis train` excludes > 200 kb | 17 → **2,402 loci** |

v3-active includes BGC0001436 and BGC0002116 (absent from the 4.0 JSON) and BGC0002224 (retired in 4.0), which are
active in 5.0rc1. The negatives are pool v2 (objective 01), whose isolate test backbones are held out of training.

## Test backgrounds (`03_background.sh`)

The shuffled set uses 124 genomes that are not in pool v2
(`../00-data-collection/accessions/testset_background_genomes.tsv`, selected by `workers/sample_background.py`):

| | Genomes | Selection |
|---|---|---|
| prokaryotes | 102 (85 Bacteria, 17 Archaea) | GTDB r226 species representatives of no GTDB species of a pool assembly, `balanced_sample` over phyla (seed 43) |
| eukaryotes | 22 (8 fungi, 14 protozoa) | RefSeq, pool assemblies and species removed, `sample_refseq` per phylum (seed 1), 22 drawn per group in proportion |

Eukaryote GBFFs are cleaned first (`cleaned/cleaning.tsv`): organelle records are dropped, only source + CDS features
are kept, and a CDS repeating an earlier location or name (locus_tag, else gene, else protein_id; first kept) or
with overlapping exons is dropped (antiSMASH 8.0.1 rejects such records). Processing is then `sanntis
build-negatives`' own code: pyrodigal (single; meta < 100 kb) on prokaryotes, the cleaned RefSeq CDS on eukaryotes,
antiSMASH 8.0.1 `--minimal`, native mask, BGC-free runs. InterProScan is not run on the backgrounds: the test contigs
are nucleotide GenBank and SanntiS annotates their proteins when scoring.

## Test sets (`04_testsets.sh`)

Test B = objective 01's archives with accession >= BGC0002750, active or pending, locus <= 200 kb: **375 BGCs** (281
prokaryotic, 94 eukaryotic: 87 fungi, 3 plants, 4 other eukaryotes; BGC0002977, a 4.2 Mb pending locus, is the one
candidate excluded). Test A = its 178 active entries. Kingdom and lineage from the MIBiG NCBI taxid (Eukaryota →
euk, split by the NCBI lineage; anything else → prokaryotes, incl. the giant-virus entry BGC0002846). Each BGC is
inserted once per set, on the + strand, with the CDS SanntiS reads from the archive (`/pseudo` dropped):

| Set | Contigs | Backbones | Rule (seed 0 + crc32(track)) |
|---|---|---|---|
| isolate | 70 (55 proka + 15 euk), 5-7 BGCs each | pool v2's isolate test backbones (longest BGC-free run >= 50 CDS of a genome, held out of training), sliced first CDS start → last CDS end | BGCs dealt round-robin (by accession) over the kingdom's backbones; random intergenic insertion points (the two ends and every gap between CDS reaches) |
| shuffled | 375 (281 proka + 94 euk), 1 BGC each | background contigs with >= 60 gene loci, cropped to <= 400 loci | BGCs in a random order dealt round-robin over genomes, then over each genome's contigs (random order, cycled); loci permuted with their upstream spacer, BGC inserted at a locus boundary in the middle 40-60 %; backbone CDS in an antiSMASH region flagged `native_bgc` |

Per contig: `<track>_contig_<k>.gbk` (record `synth_iso_` / `synth_shuf_<track>_<k>`; CDS with `locus_tag` =
`protein_id` = CDS id, `/translation`, lifted locations incl. eukaryotic joins), `.cds.parquet` (the truth table
objective 03 reads: `source_kind`, `entity_id`, `native_bgc`, ...), `.ground_truth.bed`. `test_bgcs.parquet` carries
Test A / B membership, kingdom, lineage and classes; SanntiS tokens and novelty are added at scoring (05).

## Scoring notes (`05_score.sh`)

`sv1.metrics` adds the test-set extras (the native-ORF rule on shuffled backgrounds, the masked bootstrap, AUPRC) on
top of `sanntis.train.metrics`. Test A and Test B share the contigs (Test A only masks the pending ground truth), so
precision is identical between them by construction. On the shuffled set, predictions over a backbone's native
antiSMASH region and no test BGC are ignored. Run the scorer as a SLURM step (the runbook default); peak memory is
~3 GB.

## Expected (`expected/`)

`v3active_manifest.tsv`: the training archives (above). Frozen from the run: TODO (`background_counts.tsv`, the
bundle report tables and `test_*.tsv`).
