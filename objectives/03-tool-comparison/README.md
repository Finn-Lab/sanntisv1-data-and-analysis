# 03 · Tool comparison

**Objective:** compare SanntiS 1.2.0 with other BGC tools, every method on the same contigs and CDS, scored with the
SanntiS metric code (`sanntis.train.metrics`).

- **Detection** on objective 02's nucleotide test sets (isolate: the MIBiG test BGCs inserted into the BGC-free
  test backbones held out of negative pool v2; shuffled: into gene-shuffled contigs of genomes not in the pool; 445
  contigs, 375 test BGCs): antiSMASH 7.1.0 lineage-matched (bacteria on prokaryotic contigs, fungiSMASH on
  eukaryotic ones, plantiSMASH 1.0 over plant BGCs; strict / relaxed / loose), GECCO 0.9.5 and 0.11.0, SanntiS
  0.9.5.0, BGC-Prophet 0.1.2, and SanntiS 1.2.0 (the v3-active model of objective 02).
- **Product classification:** BGCat vs the SanntiS NPC classifier (NPClassifier pathway / superclass) on the test
  BGCs.

```bash
bash objectives/03-tool-comparison/01_inputs.sh           # objective 02's test sets and evaluation ($DATA_ROOT/02)
bash objectives/03-tool-comparison/02_tools.sh [run ...]  # the tools (containers), incl. ips576 for SanntiS 0.9.5
bash objectives/03-tool-comparison/03_sanntis_tokens.sh   # SanntiS tokens of the test proteins (objective 02's)
bash objectives/03-tool-comparison/04_score.sh            # detection metrics (+ comparison with expected/ once frozen)
bash objectives/03-tool-comparison/05_classification.sh   # BGCat vs SanntiS NPC
bash objectives/03-tool-comparison/06_report.sh           # report.md + figures
```

Status: TODO (filled from the run).

## Inputs

Objective 02's test sets (`$DATA_ROOT/02/testsets`, published as `nucleotide_testsets_v3`) and evaluation
(`$DATA_ROOT/02/eval`: the v3-active SanntiS model, the test-BGC table with Test A / B, classes, lineage, novelty,
tokens; the NPClassifier labels; the SanntiS tokens of the test proteins; the Test A / B scores; published as
`v3active_eval_v3`).

Also fetched from their publishers (pinned, sha256-checked; never redistributed by this repository or its deposits):
the BGC-Prophet 0.1.2 models (GitHub attachment `files/12733164/model.tar.gz`) and the ESM2 8M weights it loads
(`dl.fbaipublicfiles.com/fair-esm`), the ESM C 600M weights for BGCat (Hugging Face, commit `e4d83bc7`), and the
NPClassifier ontology (`index_v1.json`, GitHub commit `78e52f1d`).

## Steps

| Script | Does | Output (`$DATA_ROOT/03/…`) |
|---|---|---|
| `01_inputs.sh` | links objective 02's test sets and evaluation | `inputs/{nucleotide_testsets, v3active_eval}` |
| `02_tools.sh` | runs: `gecco-0.9.5`, `gecco-0.11.0`, `ips576` (InterProScan 5.76-107.0 per contig, protein ids = `/protein_id`), `sanntis-0.9.5`, `bgc-prophet`, `antismash-{relaxed,strict,loose}`, `plantismash`, `bgcat` | `tools/<run>/<set>_<track>/`, `ips576/<set>/<track>/` |
| `03_sanntis_tokens.sh` | SanntiS tokens of the test proteins: objective 02's (InterProScan 6.0.3.0 / InterPro 107.0 via SanntiS, md5-keyed) | `sanntis/tokens.parquet` |
| `04_score.sh` | metrics with 95% CIs (bootstrap over source contigs), strata, PR curves, for SanntiS 1.2.0 and every run under `tools/`; subset scopes with CIs: `proka` / `euk` (contigs of that kingdom) and, on the shuffled sets only, `fungi` (contigs whose one BGC is fungal), and `paired.tsv`: SanntiS − antiSMASH relaxed on the `euk` / `fungi` scopes, paired bootstrap over contigs; compares with `expected/` when its tables exist | `score/` |
| `05_classification.sh` | SanntiS NPC vs BGCat (needs `02_tools.sh bgcat`): top-1 pathway, top-1 / top-3 superclass, coverage | `classification/` |
| `06_report.sh` | report | `report/report.md`, `report/figures/` |

How each tool is run (`02_tools.sh`): every tool gets the contig GenBank. antiSMASH 7.1.0 `--genefinding-tool none
--hmmdetection-strictness <level>`, `--taxon bacteria` / `fungi` by track; plantiSMASH 1.0 `--taxon plants
--genefinding none` on eukaryotic contigs; GECCO `gecco run -g` (own gene calls, defaults); SanntiS 0.9.5.0
`--ip-file <IPS 5.76 TSV>`; BGC-Prophet on each contig's proteins in gene order (`--threshold 0.5 --max_gap 3
--min_count 2`, the README values); BGCat on the proteins of the antiSMASH 7.1.0 (relaxed) `biosynthetic` /
`biosynthetic-additional` genes inside each test BGC, as its README prescribes for antiSMASH users.

## Containers (`env/containers.tsv`)

GECCO 0.9.5 / 0.11.0 and SanntiS 0.9.5.0 run from digest-pinned biocontainers. Four images are ours, defined in
`env/docker/` (never pushed to a registry: no third-party code or weights are redistributed; point `IMAGE_<tool>` at
a local image, or `IMAGE_<tool>=host` with `<TOOL>_ENV`, as below):

| Image | Why |
|---|---|
| `antismash-7.1.0` | the biocontainer ships no databases and is read-only; the layer downloads + prepares the 7.1.0 database set |
| `plantismash-1.0` | no container exists (Python 2.7, antiSMASH 3) |
| `bgc-prophet-0.1.2` | the biocontainer cannot run the pipeline: fair-esm 1.0.2 has no ESM2, and Python 3.14's forkserver cannot pickle BGC-Prophet's LMDB dataset |
| `bgcat` | no container; `embed.py` hard-codes a second GPU (`cuda:1`), patched to the first GPU or the CPU |

## Tool runs (Codon)

The runs use the following images and environments (`IMAGE_<tool>` / `<TOOL>_ENV` / `<TOOL>_HOME`, see
`02_tools.sh`):

| Run | Version / tag | Image or environment |
|---|---|---|
| `gecco-0.9.5`, `gecco-0.11.0`, `sanntis-0.9.5` | 0.9.5, 0.11.0, 0.9.5.0 (reports 0.9.4.1) | the pinned biocontainers (digests in `env/containers.tsv`), pulled once to SIF |
| `ips576` | InterProScan 5.76-107.0 | `quay.io/microbiome-informatics/interproscan:5.76-107.0` + the 5.76-107.0 data release mounted at `/opt/interproscan/data` |
| `antismash-{relaxed,strict,loose}` | 7.1.0 | `quay.io/biocontainers/antismash:7.1.0` as a SIF, databases of `download-databases` 7.1.x (Pfam 35.0) |
| `plantismash` | 1.0 (tag `1.0` = `869a16b`) | host conda env (Python 2.7) + checkout |
| `bgc-prophet` | 0.1.2 | host venv as its Dockerfile: Python 3.12.13, bgc-prophet 0.1.2, fair-esm 2.0.0, torch 2.5.1+cpu |
| `bgcat` | `4f77517` | host venv as its Dockerfile: Python 3.12.13, esm 3.2.3, torch 2.15.0, scikit-learn 1.9.1, rdkit 2026.9.2; checkout with the device patch; ESM C 600M `e4d83bc7`; GPU |

GECCO 0.9.5 and 0.11.0 embed the same CRF model (`gecco/crf/model.pkl.md5` `e7066c9d…`), so their regions are
identical. BGC-Prophet keeps one region per 128-gene window, which costs it recall
on the isolate sets (several BGCs per contig).

## Expected (`expected/`)

Frozen from the run: TODO (`summary.tsv`, `metrics.tsv` with the proka / euk / fungi scopes and CIs, `strata.tsv`,
`radar.tsv`, `curve.tsv`, `auprc.tsv`, `paired.tsv`, `paired_bgc.tsv`, `classification/`); `04_score.sh` and
`05_classification.sh` compare every one.

## Caveats for the paper

- **BGCat's training data may include the test BGCs.** Its paper trains on MIBiG 3.1 for the temporal split, but
  does not say what the released weights were trained on. Report the rows of test BGCs that are *not in MIBiG 4.0*.
- **BGC-Prophet keeps one region per 128-gene window** (its `output` step keeps the longest run of BGC genes), so on
  the isolate contigs, which carry several BGCs, its recall is low while on the shuffled contigs (one BGC each) it
  is close to the others. That is the tool's documented behaviour with its README parameters.
