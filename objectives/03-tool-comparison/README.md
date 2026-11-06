# 03 · Tool comparison

**Objective:** compare SanntiS 1.0 with other BGC tools, every method on the same contigs and CDS, scored with the
SanntiS 1.0 metric code (`sanntis.train.metrics`).

- **Detection** on objective 02's nucleotide test sets (isolate: the MIBiG test BGCs inserted into BGC-free stretches
  of genomes disjoint from the training negatives; shuffled: into gene-shuffled contigs of such genomes; 445
  contigs, 375 test BGCs): antiSMASH 7.1.0 lineage-matched (bacteria on
  prokaryotic contigs, fungiSMASH on eukaryotic ones, plantiSMASH 1.0 over plant BGCs; strict / relaxed / loose),
  GECCO 0.9.5 and 0.11.0, SanntiS 0.9.5.0, BGC-Prophet 0.1.2, and SanntiS 1.0 (the v3-active model of objective 02).
- **Product classification:** BGCat vs the SanntiS 1.0 NPC classifier (NPClassifier pathway / superclass) on the
  test BGCs.

```bash
bash objectives/03-tool-comparison/01_inputs.sh           # test sets + objective 02's evaluation (deposits or $DATA_ROOT/02)
bash objectives/03-tool-comparison/02_tools.sh [run ...]  # the tools (containers), incl. ips576 for SanntiS 0.9.5
bash objectives/03-tool-comparison/03_sanntis_tokens.sh   # SanntiS tokens (objective 02's, or InterProScan 6 via SanntiS)
bash objectives/03-tool-comparison/04_score.sh            # detection metrics (+ comparison with expected/ once frozen)
bash objectives/03-tool-comparison/05_classification.sh   # BGCat vs SanntiS NPC
bash objectives/03-tool-comparison/06_report.sh           # report.md + figures
```

Status (2026-11-04): on objective 02's v2 test sets (backgrounds disjoint from pool v1): SanntiS 1.0 (v3-active,
t* 0.25) and every tool re-run on them (2026-11-02, "Tool runs" below); `expected/` frozen from that run.

## Inputs (`env/zenodo.tsv`)

| Deposit | Content |
|---|---|
| `nucleotide_testsets_v2` | objective 02's test sets, per contig: GenBank with CDS calls (the input of every tool), CDS truth table, BGC BED |
| `v3active_eval_v2` | objective 02's evaluation: the v3-active SanntiS model, the test-BGC table (Test A / B, classes, lineage, novelty, tokens), the NPClassifier labels, the SanntiS tokens of the test proteins, the Test A / B scores. With an objective 02 run in `$DATA_ROOT` (or `V3ACTIVE_EVAL=<its eval dir>`) that is used instead |

Also fetched from their publishers (pinned, sha256-checked; never redistributed by this repository or its deposits):
the BGC-Prophet 0.1.2 models (GitHub attachment `files/12733164/model.tar.gz`) and the ESM2 8M weights it loads
(`dl.fbaipublicfiles.com/fair-esm`), the ESM C 600M weights for BGCat (Hugging Face, commit `e4d83bc7`), and the
NPClassifier ontology (`index_v1.json`, GitHub commit `78e52f1d`).

## Steps

| Script | Does | Output (`$DATA_ROOT/03/…`) |
|---|---|---|
| `01_inputs.sh` | the test sets (deposit) and objective 02's evaluation (`$DATA_ROOT/02/eval`, `V3ACTIVE_EVAL`, or the deposit) | `inputs/{nucleotide_testsets, v3active_eval}` |
| `02_tools.sh` | runs: `gecco-0.9.5`, `gecco-0.11.0`, `ips576` (InterProScan 5.76-107.0 per contig, protein ids = `/protein_id`), `sanntis-0.9.5`, `bgc-prophet`, `antismash-{relaxed,strict,loose}`, `plantismash`, `bgcat` | `tools/<run>/<set>_<track>/`, `ips576/<set>/<track>/` |
| `03_sanntis_tokens.sh` | SanntiS tokens of the 96,278 test proteins: copied from `v3active_eval/tokens.parquet` (md5-keyed), else IPS 6.0.3.0 / InterPro 107.0 via SanntiS | `sanntis/tokens.parquet` |
| `04_score.sh` | metrics with 95% CIs (bootstrap over source contigs), strata, PR curves, for SanntiS 1.0 and every run under `tools/`; subset scopes with CIs: `proka` / `euk` (contigs of that kingdom) and, on the shuffled sets only, `fungi` (contigs whose one BGC is fungal), and `paired.tsv`: SanntiS 1.0 − antiSMASH relaxed on the `euk` / `fungi` scopes, paired bootstrap over contigs ([F-03-4](../../docs/FINDINGS.md#f-03-4-eukaryotic-and-fungal-bgcs-precision-by-attributable-scope-verified-2026-11-03)); compares with `expected/` when its tables exist | `score/` |
| `05_classification.sh` | SanntiS NPC vs BGCat (needs `02_tools.sh bgcat` on the v2 sets): top-1 pathway, top-1 / top-3 superclass, coverage | `classification/` |
| `06_report.sh` | report | `report/report.md`, `report/figures/` |

How each tool is run (`02_tools.sh`): every tool gets the contig GenBank. antiSMASH 7.1.0 `--genefinding-tool none
--hmmdetection-strictness <level>`, `--taxon bacteria` / `fungi` by track; plantiSMASH 1.0 `--taxon plants
--genefinding none` on eukaryotic contigs; GECCO `gecco run -g` (own gene calls, defaults); SanntiS 0.9.5.0
`--ip-file <IPS 5.76 TSV>`; BGC-Prophet on each contig's proteins in gene order (`--threshold 0.5 --max_gap 3
--min_count 2`, the README values); BGCat on the proteins of the antiSMASH 7.1.0 (relaxed) `biosynthetic` /
`biosynthetic-additional` genes inside each test BGC, as its README prescribes for antiSMASH users.

## Containers (`env/containers.tsv`)

GECCO 0.9.5 / 0.11.0 and SanntiS 0.9.5.0 run from digest-pinned biocontainers. Four images are ours, defined in
`env/docker/` (built locally only, never pushed to a registry: no third-party code or weights are redistributed;
point `IMAGE_<tool>` at a local image, or `IMAGE_<tool>=host` with `<TOOL>_ENV`, as the v2 runs did: "Tool runs"
below):

| Image | Why |
|---|---|
| `antismash-7.1.0` | the biocontainer ships no databases and is read-only; the layer downloads + prepares the 7.1.0 database set |
| `plantismash-1.0` | no container exists (Python 2.7, antiSMASH 3) |
| `bgc-prophet-0.1.2` | the biocontainer cannot run the pipeline: fair-esm 1.0.2 has no ESM2, and Python 3.14's forkserver cannot pickle BGC-Prophet's LMDB dataset |
| `bgcat` | no container; `embed.py` hard-codes a second GPU (`cuda:1`), patched to the first GPU or the CPU |

## Tool runs on the v2 sets (Codon)

The `env/docker/` images were not built; the runs below use existing images and environments, recorded here
(`IMAGE_<tool>` / `<TOOL>_ENV` / `<TOOL>_HOME`, see `02_tools.sh`). Dates are run starts (UTC).

| Run | Version / tag | Image or environment | Run |
|---|---|---|---|
| `gecco-0.9.5` | 0.9.5 | `quay.io/biocontainers/gecco:0.9.5--pyhdfd78af_0` (digest in `env/containers.tsv`) | 2026-11-02 |
| `gecco-0.11.0` | 0.11.0 | `quay.io/biocontainers/gecco:0.11.0--pyhdfd78af_0` (digest in `env/containers.tsv`) | 2026-11-02 |
| `ips576` | InterProScan 5.76-107.0 | `quay.io/microbiome-informatics/interproscan:5.76-107.0` (created 2025-10-27) + the 5.76-107.0 data release mounted at `/opt/interproscan/data` | 2026-11-02 (+ helper array, below) |
| `sanntis-0.9.5` | 0.9.5.0 (reports 0.9.4.1) | `quay.io/biocontainers/sanntis:0.9.5.0--pyhdfd78af_0` (digest in `env/containers.tsv`) | 2026-11-02 |
| `antismash-{relaxed,strict,loose}` | 7.1.0 | `quay.io/biocontainers/antismash:7.1.0` converted to a writable sandbox / SIF (2026-10-19), databases of `download-databases` 7.1.x (Pfam 35.0; 2023-11-24): the image and databases of the archived runs | 2026-11-02 |
| `plantismash` | 1.0 (tag `1.0` = `869a16b`) | host conda env (Python 2.7) + checkout | 2026-11-02 |
| `bgc-prophet` | 0.1.2 | host venv: Python 3.12.13, bgc-prophet 0.1.2, fair-esm 2.0.0, torch 2.5.1+cpu, lmdb 3.0.0 (as its Dockerfile) | 2026-11-02 |
| `bgcat` | `4f77517` | host venv: Python 3.12.13, esm 3.2.3, torch 2.14.1, rdkit 2026.9.1; checkout with the Dockerfile's device patch; ESM C 600M `e4d83bc7` (sha256 `8ef856e1…`) | 2026-11-02 |

`ips576` ran as one 32-core step (8 contigs × 4 cores) plus, from 10:57 UTC, a 24-task SLURM array
(`sbatch --array=0-23 -c 4`) running the same `ips576_one` on contigs taken from the other end of the list, each claimed
with a `mkdir` lock and written under a per-process temporary name; the step skips contigs that already have an
output, and the tables do not depend on which process wrote them. 383 of the 445 came from the helpers; all were done
at 12:34 UTC. The `ips576` and first `gecco-0.9.5` step jobs ended FAILED after writing their markers and outputs:
`02_tools.sh` was edited while they ran and bash re-reads a running script. `gecco-0.11.0` ran twice (first attempt:
3 contigs failed on parallel `docker://` pulls; re-run from a pre-pulled SIF). SanntiS 0.9.5 ran 12:35–12:47 UTC.
GECCO 0.9.5 and 0.11.0 call identical regions on all 369 contigs with predictions (their rows in the report are
the same): both embed the same CRF model (`gecco/crf/model.pkl.md5` `e7066c9d…`); only the output columns (`bgc_id` →
`cluster_id`) and float rounding differ.
BGC-Prophet keeps one region per 128-gene window, which costs it recall on the isolate sets (several BGCs per
contig): see [F-03-1](../../docs/FINDINGS.md#f-03-1-bgc-prophet-reports-at-most-one-bgc-per-128-gene-window-verified-2026-11-02).

InterProScan 5.76: an earlier host install on Codon is gone. On a 98-protein contig, this image + data gives rows
identical to that install's for 12 of the 14 member databases (Pfam included); PIRSF E-values and a few Gene3D domain
boundaries differ, both dependent on how many proteins share a run (that install: 10,000-protein chunks; here one run
per contig). See [F-03-3](../../docs/FINDINGS.md).

## Expected (`expected/`)

Frozen 2026-11-04 from the `04_score.sh` run of 2026-11-03 (`expected/PROVENANCE.json`, `expected/NOTE.txt`):
`summary.tsv`, `metrics.tsv` (with the proka / euk / fungi scopes and CIs), `strata.tsv`, `radar.tsv`, `curve.tsv`,
`auprc.tsv`, `paired.tsv`; `04_score.sh` compares every one. A second, independent `score.py` run gave byte-identical
tables. `expected/classification/`: the `05_classification.sh` tables (frozen 2026-11-04, compared by that step; they
hold for re-runs from the same BGCat predictions, which ran on a GPU). `expected/v1sets/`: the tables of the previous
sets from the pre-split release run (t* 0.35; `summary.tsv`, `metrics.tsv` with CIs and per kingdom, `strata.tsv`,
`radar.tsv`, `curve.tsv`, `auprc.tsv`: SanntiS 1.0, antiSMASH 7.1.0 and GECCO 0.9.5 rows, e.g. A-isolate F1_macro_T
0.802 / relaxed 0.624 / 0.645), the ones the verification below reproduced.

## Verification on the previous sets (Codon, 2026-10-23)

- `04_score.sh` from the deposits: every expected table identical (1,680 metric rows incl. CIs, 250 radar,
  28 strata, 72 curve, 16 AUPRC, 20 summary rows).
- `02_tools.sh gecco-0.9.5` (pinned biocontainer): all 445 contigs, region coordinates and scores identical to the
  archived runs (828 regions).
- `02_tools.sh antismash-relaxed` with the image and databases of the archived runs (`IMAGE_antismash`,
  `ANTISMASH_DB`), 2 contigs per set and track: regions identical. The `env/docker/antismash-7.1.0` image itself
  is not built yet.
- GECCO 0.11.0 calls the same regions as 0.9.5 on these sets (scores differ in the last digit, gene numbering).
- SanntiS 0.9.5.0: the bioconda package `sanntis-0.9.5.0` reports version 0.9.4.1 (its version string was not
  bumped).
- BGC-Prophet and BGCat ran in Python environments equal to their Dockerfiles (`IMAGE_<tool>=host`); BGCat on a GPU.
- End to end from the deposits only (`01` → `03` → `04` → `05` → `06`, fresh `DATA_ROOT`): exit 0, expected tables
  identical, every new tool scored, report + 12 figures written.
- `05_classification.sh`: the SanntiS NPC rows equal objective 02's expected true-loci rows (Test B: n 322,
  pathway top-1 0.879, superclass top-1 0.738, top-3 0.867).

## Caveats for the paper

- **BGCat's training data may include the test BGCs.** Its paper trains on MIBiG 3.1 for the temporal split, but
  does not say what the released weights were trained on. On the test BGCs it covers BGCat leads on superclass
  top-1 (0.815 vs 0.752); on the 123 covered BGCs that are *not in MIBiG 4.0* it drops to 0.664 (SanntiS 0.776),
  pathway top-1 0.797 vs 0.911. Report the "not in MIBiG 4.0" rows, or ask the authors.
- **BGC-Prophet keeps one region per 128-gene window** (its `output` step keeps the longest run of BGC genes), so on
  the isolate contigs, which carry several BGCs, its recall is low (A-isolate BGC recall 0.47) while on the shuffled
  contigs (one BGC each) it is close to the others (B-shuffled F1_macro_T 0.777). That is the tool's documented
  behaviour with its README parameters, not a parsing artefact.
