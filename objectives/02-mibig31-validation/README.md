# 02 · MIBiG 3.1 validation

**Objective:** train on v3-active (MIBiG 3.1-era entries, accession < BGC0002750, active / pending in 5.0rc1: 2,419
BGCs) and score the model on Test A (5.0rc1 active, new since v3: 178) and Test B (active + pending: 375), each
inserted into isolate and gene-shuffled genomic backgrounds **whose genomes are not in the training negative pool**;
strata by class, lineage and novelty.

```bash
bash objectives/02-mibig31-validation/01_v3active.sh     # v3-active MIBiG-format archives (from objective 01's)
bash objectives/02-mibig31-validation/02_train.sh        # sanntis train → the v3-active model
bash objectives/02-mibig31-validation/03_background.sh   # test-background pool (gene calls, antiSMASH 8.0.1, strip)
bash objectives/02-mibig31-validation/04_testsets.sh     # isolate + shuffled nucleotide test sets, Zenodo deposit
bash objectives/02-mibig31-validation/05_score.sh        # IPS6 tokens of the test proteins, then Test A / B scoring
bash objectives/02-mibig31-validation/06_report.sh       # report.md + figures, the v3active_eval_v2 deposit
```

Status: all six steps run on Codon (2026-10-29 / 30): training 40 min (head job + SLURM fan-out), tokens 2.5 h,
scoring 4 min. Headline (v3-active model, t* 0.25): F1_macro_T A-isolate 0.754, B-isolate 0.765, A-shuffled 0.727,
B-shuffled 0.743; BGC precision 0.960 (isolate) / 0.933 (shuffled); see `report/report.md`.

## Steps

| Script | Does | Output (`$DATA_ROOT/02/…`) |
|---|---|---|
| `01_v3active.sh` | worker `v3_archives.py`: the < BGC0002750 members of objective 01's archives, verbatim; checked with SanntiS's external-data reader and against `expected/v3active_manifest.tsv` | `inputs/mibig_v3-active/` |
| `02_train.sh` | `sanntis train --mibig none --external-gbk/-json` on v3-active + pool v1; checks the collected pool (2,402 loci) | `run/` (`bundle/`, `report/`, `tune/`) |
| `03_background.sh` | worker `background_pool.py`: objective 00's 124 test-background genomes through `sanntis build-negatives`' stages up to the stripping; no InterProScan | `background/background/*.parquet` |
| `04_testsets.sh` | worker `build_testsets.py`: the Test B BGCs inserted into the background pool as isolate and gene-shuffled nucleotide contigs; the deposit tarball (`env/zenodo.tsv`: `nucleotide_testsets_v2`) | `testsets/{testsets/, test_bgcs.parquet, contigs.tsv, manifest.json}`, `deposit/` |
| `05_score.sh` | step `02_tokens`, worker `tokens.py`: InterProScan 6 through SanntiS's annotation (md5 cache) → tokens + Pfam sets per test protein; step `02_score`, worker `score_testsets.py`: Test A / B on both sets at the bundle's t* with `sanntis.train.metrics` (CIs, strata, radar, PR curves, NPC), novelty vs the training positives, the test BGCs' NPC labels from objective 01's `collect`; `expected/release_v1sets/` reported side by side (`score/reference.tsv`), `expected/` compared | `eval/{tokens.parquet, test_bgcs.parquet, labels/, bundle → run/bundle, score/{A,B}-{isolate,shuffled}.pkl, score/test_*.tsv}` (`V3ACTIVE_EVAL` of objective 03) |
| `06_report.sh` | step `02_report`, worker `report_v3active.py`: training and OOF-validation tables, Test A / B metrics with CIs, threshold sweeps (reference sweep dashed), strata radars, novelty, NPC; `expected/release_v1sets/` as reference columns, `expected/` compared; step `02_deposit`: the `v3active_eval_v2` deposit (bundle, test-BGC table, labels, tokens, scores) objective 03 reads | `report/{report.md, summary.tsv, figures/}`, `deposit/sanntisv1_v3active_eval_v2.tar.gz` |

## Training data (`01_v3active.sh`)

| | BGCs |
|---|---|
| objective 01 archives (5.0rc1 active + pending) with accession < BGC0002750 | 2,419 (2,418 active, 1 pending) |
| of which loci from nuccore | 18 (incl. the 15 at their changed 5.0rc1 coordinates) |
| `sanntis train` excludes > 200 kb | 17 → **2,402 loci** |

The previous v3-active (2,416) also required the entry to be active in MIBiG 4.0; BGC0001436 and BGC0002116 (absent
from the 4.0 JSON) and BGC0002224 (retired in 4.0) are active in 5.0rc1 and are included here (decision 2026-10-29).

## Test backgrounds (`03_background.sh`)

The previous test sets (objective 03's v5 isolate / v5shuf2, and the `expected/release_v1sets/` tables) used pool v1 genomes as
backgrounds: the isolate backbones were reserved runs of pool genomes and the shuffled contigs came from every pool
genome, i.e. the training negatives. Here the backgrounds are 124 other genomes
(`../00-data-collection/accessions/testset_background_genomes.tsv`, selected by `workers/sample_background.py`):

| | Genomes | Selection |
|---|---|---|
| prokaryotes | 102 (85 Bacteria, 17 Archaea; 102 phyla) | GTDB r226 species representatives not in pool v1, `balanced_sample` over phyla (seed 43) |
| eukaryotes | 22 (9 fungi, 13 protozoa) | RefSeq, pool v1 assemblies and species removed, `sample_refseq` per phylum (seed 1), 22 drawn per group in proportion |

No assembly or species is shared with pool v1 (11 of the 102 GTDB genera also occur in it). Eukaryote GBFFs are
cleaned first (`cleaned/cleaning.tsv`): antiSMASH 8.0.1 rejected 6 of the 22 RefSeq annotations as given (CDS
repeating another's location or locus_tag from alternative transcripts, CDS with overlapping exons, duplicate gene
names on a mitochondrial record), so organelle records are dropped, only source + CDS features are kept, and a CDS
repeating an earlier location or name (locus_tag, else gene, else protein_id; first kept: one CDS per gene, as
objective 01 does for nuccore loci) or with overlapping exons is dropped.
Processing is then `sanntis build-negatives`' own code: pyrodigal (single; meta < 100 kb) on prokaryotes, the
cleaned RefSeq CDS on eukaryotes,
antiSMASH 8.0.1 `--minimal`, native mask, BGC-free runs and the isolate backbones (55 prokaryotic + 15 eukaryotic,
longest run >= 50 CDS, the pool's §2.5 rule). InterProScan is not run on the pool: the test contigs are nucleotide
GenBank and SanntiS annotates their proteins when scoring.

## Test sets (`04_testsets.sh`)

Test B = objective 01's archives with accession >= BGC0002750, active or pending, locus <= 200 kb: **375 BGCs** (281
prokaryotic, 94 eukaryotic: 87 fungi, 3 plants, 4 other eukaryotes; BGC0002977, a 4.2 Mb pending locus, is the one
candidate excluded). Test A = its 178 active entries. Kingdom and lineage follow the release run (Eukaryota →
euk, split by the NCBI lineage; anything else → prokaryotes, incl. the giant-virus entry BGC0002846). Each BGC is
inserted once per set, on the + strand, with the CDS SanntiS reads from the archive (`/pseudo` dropped):

| Set | Contigs | Backbones | Rule (sanntis_dev builders / release feature-level sets, seed 0 + crc32(track)) |
|---|---|---|---|
| isolate | 70 (55 proka + 15 euk), 5-7 BGCs each | the pool's test backbones (longest BGC-free run >= 50 CDS of a genome), sliced first CDS start → last CDS end | BGCs dealt round-robin (by accession) over the kingdom's backbones; random intergenic insertion points (dev rule: the two ends and every gap between CDS reaches) |
| shuffled | 375 (281 proka + 94 euk), 1 BGC each | contigs with >= 60 gene loci of 62 / 102 prokaryotic and 21 / 22 eukaryotic genomes (the others are too fragmented), cropped to <= 400 loci | BGCs in a random order dealt round-robin over genomes, then over each genome's contigs (random order, cycled); loci permuted with their upstream spacer, BGC inserted at a locus boundary in the middle 40-60 %; backbone CDS in an antiSMASH region flagged `native_bgc` (2,234 CDS) |

Per contig: `<track>_contig_<k>.gbk` (record `synth_iso_` / `synth_shuf_<track>_<k>`; CDS with `locus_tag` =
`protein_id` = CDS id, `/translation`, lifted locations incl. eukaryotic joins), `.cds.parquet` (the truth table
objective 03 reads: `source_kind`, `entity_id`, `native_bgc`, ...), `.ground_truth.bed`. `test_bgcs.parquet` carries
Test A / B membership, kingdom, lineage and classes; SanntiS tokens and novelty are added at scoring (05).
Every lifted CDS re-translates to its qualifier (44,580 of 44,743 in a trial; the rest are RefSeq 3'-partial
CDS with one extra residue and ambiguous codons pyrodigal wrote as X). 460 MB on disk, 164 MB deposit.

## Scoring notes (`05_score.sh`)

`sv1.metrics` restores the test-set extras SanntiS dropped before 1.0.0b1 (the native-ORF rule on shuffled
backgrounds, the masked bootstrap, AUPRC) on top of `sanntis.train.metrics`. Test A and Test B share the contigs
(Test A only masks the pending ground truth), so precision is identical between them by construction. On the
shuffled set, predictions over a backbone's native antiSMASH region and no test BGC are ignored (56 at t*); counting
them as background gives Precision_T 0.566 / BGC precision 0.811 instead of 0.658 / 0.933. Run the scorer as a SLURM
step (the runbook default): the login node killed a direct run; peak memory is ~3 GB.

## Expected (`expected/`)
`v3active_manifest.tsv`: the training archives (above); `background_counts.tsv`: the test backgrounds. The other
tables are frozen from this objective's run (2026-11-02, `PROVENANCE.json`: SanntiS v1.0.0-beta.1, 2,402 loci,
t* 0.25): `tune.tsv`, `validation_oof.tsv` (from `tune/tune.json`), `validation_strata.tsv`,
`validation_class_curves.tsv`, `npc_validation.tsv`, `sweep.tsv` (the bundle's report) and `test_*.tsv` (Test A / B
on the v2 sets, e.g. F1_macro_T A-isolate 0.754, A-shuffled 0.727, B-isolate 0.765, B-shuffled 0.743).
`05_score.sh` compares its `score/test_*.tsv` and `06_report.sh` the bundle's report tables with them
(`sv1.expected compare`).

`release_v1sets/`: the same tables from the pre-split release run (v3-active of 2,416 BGCs, t* 0.35, test sets on
pool v1 backgrounds: F1_macro_T A-isolate 0.795, A-shuffled 0.785, B-isolate 0.818, B-shuffled 0.811), the
reference columns and dashed sweeps of the report; a reference, not a target.
