# Findings

One place for the insights, caveats and decisions that come out of the runs, written so they can be quoted in reports
and reused elsewhere. Each entry is a claim with its evidence (run, files, numbers) and its scope: what it supports and
what it does **not**. Objective READMEs keep the how-to; this file keeps what the results mean.

Conventions: newest first within a section; IDs (`F-<objective>-<n>`) are stable, so cite them; status is
`verified` (re-checked from the outputs), `open` (needs more work) or `decision` (a call made on a verified result).

## Objective 03: tool comparison

### F-03-5 BGC-level headline: SanntiS 1.0 vs every tool, paired (verified 2026-11-05)

**Claim.** On B-shuffled (375 BGCs), SanntiS 1.0 finds significantly more test BGCs than every other tool, and its
BGC-level F1 (`f1_macro`: BGC precision with class-averaged BGC recall) is significantly higher than SanntiS 0.9.5,
BGC-Prophet and GECCO, but **not** than antiSMASH 7.1.0 relaxed (lineage-matched), whose BGC precision is higher.

**Evidence** (`objectives/03-tool-comparison/workers/paired_bgc.py`, step `03_paired_bgc` of `04_score.sh`, output
`data/03/score/paired_bgc.tsv`, frozen in `expected/`; same inputs, regions and `score.paired` as `score.py`, 1,000
replicates over contigs; SanntiS 1.0 − tool, 95% CI):

| Scope | vs | BGC F1 | BGC recall | BGC precision |
|---|---|---|---|---|
| all (375) | antiSMASH 7.1.0 | +0.026 (−0.009 to +0.063), p 0.14 | +0.075 (+0.043 to +0.107) | −0.064 (−0.095 to −0.037) |
| all | SanntiS 0.9.5 | +0.083 (+0.050 to +0.120) | +0.077 (+0.051 to +0.107) | +0.048 (+0.004 to +0.091) |
| all | BGC-Prophet 0.1.2 | +0.080 (+0.040 to +0.123) | +0.117 (+0.075 to +0.160) | +0.049 (−0.003 to +0.109), p 0.08 |
| all | GECCO 0.11.0 | +0.317 (+0.255 to +0.382) | +0.144 (+0.104 to +0.184) | +0.417 (+0.329 to +0.502) |
| proka (281) | antiSMASH 7.1.0 | +0.008 (−0.026 to +0.053), p 0.67 | +0.057 (+0.025 to +0.096) | −0.069 (−0.102 to −0.037) |
| euk (94) | antiSMASH 7.1.0 | +0.080 (+0.003 to +0.187), p 0.04 | +0.128 (+0.064 to +0.202) | −0.050 (−0.091 to −0.010) |
| fungi (87) | fungiSMASH | +0.074 (−0.008 to +0.184), p 0.11 | +0.115 (+0.057 to +0.184) | −0.053 (−0.101 to −0.011) |

Gene-level F1 (`f1_macro_t`) vs BGC-Prophet on all: +0.001 (−0.038 to +0.041), p 0.93 (the gene-level tie).
The fungi / antiSMASH rows reproduce `score/paired.tsv` exactly.

**Supports.** Headline (D3.2): class-averaged BGC recall (`bgc_macro_recall`), all scope, SanntiS 1.0 0.880 vs
antiSMASH +0.093 (+0.043 to +0.152), BGC-Prophet +0.105 (+0.052 to +0.159), SanntiS 0.9.5 +0.112 (+0.067 to +0.167),
GECCO +0.194 (+0.130 to +0.265); fungi vs fungiSMASH 0.894 vs 0.731, +0.163 (+0.039 to +0.345). "Finds the most BGCs"
(all scopes, every tool); "highest BGC-level F1, significantly above its previous
release, BGC-Prophet and GECCO, level with antiSMASH"; "on fungal BGCs finds 98% vs 86% for fungiSMASH".

**Does not support.** "Outperforms antiSMASH / fungiSMASH at BGC level" (F1 not resolved; their BGC precision is
significantly higher, and inflated because the backgrounds were stripped with antiSMASH 8.0.1).

### F-03-4 Eukaryotic and fungal BGCs: precision by attributable scope (verified 2026-11-03)

**Correction first.** Until 2026-11-03 the `fungi` scope of `metrics.tsv` / `summary.tsv` (`sv1.tools.subset_metrics`)
took recall over fungal BGCs but precision over **every** prediction on eukaryotic contigs, so its precision equalled
the `euk` one; subset rows had no CIs. `score.py` now computes the subsets itself (`scopes()`): `proka` / `euk` with
bootstrap CIs over the subset's contigs, and `fungi` **only on the shuffled sets**, as the contigs whose one test BGC
is fungal. The isolate sets get no `fungi` scope: their 15 eukaryotic contigs carry 6.3 BGCs each, of mixed lineages,
so a prediction cannot be attributed to a lineage. All other tables (`all`, `proka`, `euk` values, curves, radar,
strata) are unchanged (compared with the 2026-11-02 run).

**Evidence** (`data/03/score/metrics.tsv`, report section "Eukaryotic and fungal BGCs"; Test B, 95% CIs over contigs):

| Method | B-shuffled fungal: Precision_T | Recall_T | F1_macro_T | BGC precision | B-isolate euk: Precision_T | F1_macro_T |
|---|---|---|---|---|---|---|
| SanntiS 1.0 (t* 0.25) | 0.546 (0.487–0.604) | 0.959 (0.918–0.990) | 0.675 (0.603–0.737) | 0.934 | 0.545 (0.471–0.617) | 0.700 (0.635–0.758) |
| antiSMASH 7.1.0 relaxed (fungiSMASH) | 0.444 (0.380–0.504) | 0.819 (0.739–0.893) | 0.544 (0.481–0.608) | 0.987 | 0.482 (0.413–0.547) | 0.597 (0.503–0.667) |
| GECCO 0.11.0 | 0.145 (0.109–0.190) | 0.836 (0.770–0.900) | 0.244 (0.191–0.306) | 0.233 | 0.151 (0.090–0.289) | 0.250 (0.161–0.408) |
| SanntiS 0.9.5.0 | 0.208 (0.165–0.252) | 0.905 (0.849–0.960) | 0.332 (0.275–0.387) | 0.721 | 0.281 (0.218–0.346) | 0.431 (0.352–0.503) |
| BGC-Prophet 0.1.2 | 0.428 (0.344–0.535) | 0.778 (0.695–0.866) | 0.545 (0.450–0.638) | 0.654 | 0.495 (0.387–0.649) | 0.552 (0.430–0.626) |

B-shuffled fungal: 87 contigs / BGCs; B-isolate euk: 15 contigs, 94 BGCs. Prokaryotic contigs, B-isolate Precision_T:
SanntiS 1.0 0.73, antiSMASH 0.55, GECCO 0.84, SanntiS 0.9.5 0.76, BGC-Prophet 0.90. Test A (27 fungal shuffled
contigs) gives the same ordering with wider CIs.

**Supports.**
- "On fungal BGCs SanntiS 1.0 beats fungiSMASH on Precision_T, Recall_T and F1_macro_T." Paired bootstrap
  (`score/paired.tsv`: both methods scored on the same contig resample, 1,000 replicates), SanntiS 1.0 − antiSMASH
  relaxed:

  | Scope | Precision_T | Recall_T | F1_macro_T | BGC precision |
  |---|---|---|---|---|
  | B-shuffled fungal (87 contigs) | +0.102 (+0.046 to +0.154), p < 0.001 | +0.140 (+0.071 to +0.217), p < 0.001 | +0.131 (+0.069 to +0.205), p < 0.001 | −0.053 (−0.101 to −0.011), p 0.010 |
  | A-shuffled fungal (27 contigs) | +0.167 (+0.071 to +0.264), p < 0.001 | +0.194 (+0.055 to +0.338), p < 0.001 | +0.183 (+0.097 to +0.285), p < 0.001 | −0.037 (−0.107 to 0.000), p 0.72 |
  | B-isolate eukaryotic (15 contigs) | +0.063 (−0.041 to +0.152), p 0.23 | +0.166 (+0.091 to +0.262), p < 0.001 | +0.103 (+0.000 to +0.210), p 0.05 | −0.057 (−0.143 to 0.000), p 0.24 |

  The precision gain holds on the fungal shuffled scope (both tests). On the 15 isolate eukaryotic contigs only the
  recall gain is resolved; precision is not, and F1_macro_T is borderline. **fungiSMASH has the higher region-level
  precision** (BGC precision, B-shuffled fungal): SanntiS 1.0 calls some regions that touch no BGC; fungiSMASH almost
  never does. Ties (difference exactly 0) count in both tails, which makes the A-shuffled BGC-precision p
  conservative.
- "Every method loses CDS-level precision on eukaryotic contigs, SanntiS 1.0 included (0.73 → 0.54 on B-isolate)."
  The SanntiS 1.0 drop matches its own OOF validation (F-01-1: euk Precision_T 0.565): a property of the model, not of
  these test sets.
- "Where the precision goes differs by tool": on the shuffled fungal scope (predictions touching only native backbone
  BGCs are already ignored) fungiSMASH's regions almost all hit the test BGC (BGC precision 0.99), so its lost CDS
  precision is boundary over-extension; GECCO's regions mostly hit nothing (BGC precision 0.23), so its loss is
  spurious regions. SanntiS 1.0 (0.93) is mostly extension.
- "SanntiS 1.0's main gain over 0.9.5 is on eukaryotes": fungal Precision_T 0.21 → 0.55, F1_macro_T 0.33 → 0.68 (proka
  F1_macro_T 0.71 → 0.76 on B-shuffled).

**Does not support.** Claims about plant or other eukaryotic BGCs (3 and 4 in Test B, inside the `euk` scope only);
lineage-level precision on the isolate sets; per-class conclusions with CIs. Per class (B-shuffled fungal, Recall_T, no
CIs) SanntiS 1.0 recovers terpenes at 1.00 vs fungiSMASH 0.68 (n = 25): an observation to check, not a result. A
precision advantage of SanntiS 1.0 on the isolate eukaryotic contigs (paired p 0.23). The paired test covers only
SanntiS 1.0 vs fungiSMASH; the other tools are compared by their per-method CIs.

### F-03-1 BGC-Prophet reports at most one BGC per 128-gene window (verified 2026-11-02)

**Claim.** BGC-Prophet 0.1.2's published pipeline cuts each contig into back-to-back 128-gene windows
(`bgc_prophet/command/split.py`; the last window is taken from the contig end) and its `output` step keeps only the
**longest** run of above-threshold genes per window (`output.py`, `max_range`; threshold 0.5, `--max_gap 3`,
`--min_count 2`). Any other BGC in the same window is dropped, even when its genes scored above the threshold.

**Evidence** (`data/03/tools/bgc-prophet_0.1.2/*/prophet.csv` against `data/03/inputs/nucleotide_testsets/*/*/
*.ground_truth.bed`; every "Yes" window carries exactly one gene run):

| Set | BGCs / contig | BGCs sharing a window | Hit, alone in window | Hit, sharing |
|---|---|---|---|---|
| isolate proka | 5.1 | 83% | 76% (n = 49) | 46% (n = 232) |
| isolate euk | 6.3 | 64% | 79% (n = 34) | 47% (n = 60) |
| shuffled proka | 1.0 | 0% | 81% (n = 281) | – |
| shuffled euk | 1.0 | 0% | 77% (n = 94) | – |

("Hit": any of the BGC's genes in a reported region.) Splitting at window edges matters little (1–7% of BGCs).
Precision is the same on both layouts (Test A Precision_T 0.745 isolate vs 0.743 shuffled), so the gap between the
isolate and shuffled results (Test A F1_macro_T 0.537 vs 0.725) is this output rule, not weaker detection.

**Supports.** "BGC-Prophet's output keeps one region per 128-gene window; where BGCs lie in the same window (isolate
sets) recall falls from ~0.8 to ~0.46 for the affected BGCs, while its per-window detection is comparable across
layouts."

**Does not support.** "BGC-Prophet fails on real genomes." The isolate sets pack 5–6 BGCs per contig, far denser than
real genomes on average (e.g. ~30 BGCs in ~7,500 genes of a *Streptomyces* genome); the shuffled sets (one BGC per
contig) are easier than reality. The real-world cost needs the share of BGCs with a neighbour within 128 genes in
complete annotated genomes (open, see below).

**Reporting.** Keep the published-pipeline row, labelled as such, with this caveat; use the shuffled sets as the main
comparison for BGC-Prophet.

**Open.**
- An "all runs per window" variant: the same per-gene scores and parameters, keeping every run, separates detection
  from the output rule. Needs a re-run that keeps `intermediate_prophet.npy` (~1.5 h CPU), a parser change, and
  `04_score` / `06_report` again.
- Real-world exposure: the neighbour-density count on public complete genomes.

### F-03-2 GECCO 0.9.5 and 0.11.0 give identical regions (verified 2026-11-02)

Both images (digests in `env/containers.tsv`; `gecco --version` 0.11.0 checked in the SIF) embed the same CRF model
(`gecco/crf/model.pkl.md5` = `e7066c9dcb381940e0e387ac10363625`). On all 369 contigs with predictions the region
coordinates are identical; only the output header (`bgc_id` → `cluster_id`) and float rounding differ. Their rows in
the report are therefore the same: report one GECCO row (or both, saying why they match), not two independent results.

### F-03-3 InterProScan 5.76 from the MGnify image vs an earlier host install (verified 2026-11-02)

An earlier host install of InterProScan 5.76-107.0 on Codon is gone. On a 98-protein contig, the image + 5.76-107.0
data release gives rows identical to that install's for 12 of 14 member databases (Pfam included); PIRSF E-values and
a few Gene3D boundaries differ, both dependent on how many proteins share a run (that install: 10,000-protein chunks;
here one run per contig, as for the earlier test sets). This affects only SanntiS 0.9.5's input; SanntiS 1.0 uses IPS 6.

## Objective 01: production models

### F-01-1 Precision_T just below the reference CI accepted (decision 2026-11-02)

Release run of 2026-10-30 (job 14471873), OOF validation on 2,647 BGCs: Precision_T 0.698 (CI 0.685–0.710) against
the pre-split reference 0.715 (proka 0.733 vs 0.746, reference above the CI; euk 0.565 vs 0.593, within).
F1_macro_T 0.776 vs 0.788. Recall and BGC-level values are within the CIs in every stratum. Attributed to the lower
t* (0.30 vs 0.35) on a larger pool (2,647 vs 2,297 OOF BGCs) with 8 instead of 9 epochs; it does not block the release.
Report the t* and epoch changes.

## Data deposits and licences

### F-00-1 Licences of the deposits; no third-party redistribution (decision 2026-11-04)

**Evidence.**
- The HF dataset `microbiome-informatics/sanntis-training-data` is tagged `license: cc-by-4.0` (card, 2026-11-04).
  Its `negatives/v1/genomes.parquet` (tag `negatives-v1`) holds the full GTDB r226 lineage strings
  (`d__…;p__…;…;s__…`) of the 205 prokaryotic genomes. So do the `lineage` columns of
  `objectives/00-data-collection/accessions/{pool_v1_genomes,testset_background_genomes}.tsv`.
- GTDB data is CC BY-SA 4.0 (https://gtdb.ecogenomic.org/licenses), and MIBiG is CC BY 4.0. NCBI / INSDC records
  carry no licence, and NCBI's disclaimer notes that submitters may hold IP rights in their data. InterPro / Pfam
  data is CC0.

**Decision.**
- HF dataset: CC BY 4.0 overall (MIBiG), except `negatives/`, which is CC BY-SA 4.0 because of the GTDB taxonomy.
  No data file changes, so the `negatives-v1` sha256 values stay valid. The licence text is updated on `main` only,
  and the tags never move. The HF model card stays Apache-2.0 and gains a training-data attribution section.
- In this repository, the code and our tables are Apache-2.0. The GTDB lineage columns keep CC BY-SA 4.0, and the
  MIBiG and NCBI content keeps its source terms (`NOTICE`).
- Zenodo deposits: the data deposits are CC BY 4.0 and the model deposit (`sanntis_model`, the Hub bundle) is
  Apache-2.0, as on the Hub (`licence` column of `env/zenodo.tsv`). None carries GTDB taxonomy: every text member and
  every string column of the parquet files of the data deposits (`nucleotide_testsets_v2`, `v3active_eval_v2`,
  `mibig5_train`) and of earlier, unpublished deposits was scanned for `d__` / `p__` (2026-11-04);
  `test_bgcs.parquet`'s `lineage` is only `prokaryotes` / `eukaryotes` plus `ncbi_tax_id`. Each record gets a
  `NOTICE.md` (`docs/zenodo/notices.py`, written by `docs/zenodo/stage.py`) as a file next to the tarball, not inside
  it, so the tarballs and their sha256 stay as the runbooks build and check them.
- **No third-party code or weights are redistributed.** The `third_party_weights` deposit is dropped. BGC-Prophet
  0.1.2's models and the ESM2 8M weights are downloaded by `02_tools.sh` from the publishers' URLs, pinned by sha256,
  as ESM C 600M already was. The `env/docker/` images are built locally and never pushed. Tools are referenced by
  version / tag / digest, URL and citation only.

**Does not cover.** This is the project's reading of the upstream terms, not legal advice. The NPClassifier API
answers (`npclassifier/v1/`) carry no stated data licence; they are treated as factual classifications under the
dataset's CC BY 4.0, with Kim et al. 2021 cited. If the BGC-Prophet attachment URL disappears, objective 03's
BGC-Prophet run cannot be reproduced; the sha256 is kept so that a copy from elsewhere can still be verified.

**Done 2026-11-04.** HF dataset `main` commit `6fa5f1f2` (card licence `[cc-by-4.0, cc-by-sa-4.0]`, per-folder table,
`NOTICE.md`, `LICENSES/`); HF model `main` commit `9e79818e` (attribution section); SanntiS `model_card()` writes the
same section (SanntiS commit `27d8a0e`, branch `docs/data-attribution`). Tags unchanged.

**Open.** Setting each Zenodo record's licence field to match `env/zenodo.tsv` when it is minted.
