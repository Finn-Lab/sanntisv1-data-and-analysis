# SanntiS 1.2: data and analyses

Runbooks that reproduce the data, models and analyses behind SanntiS 1.2.0. The tool itself (inference and
training) lives in [SanntiS](https://github.com/Finn-Lab/SanntiS) (release
[`v1.2.0`](https://github.com/Finn-Lab/SanntiS/releases/tag/v1.2.0)); this repository only *calls* it — the
negative pool through `sanntis plan-negatives` / `build-negatives`, training through `sanntis train`, metrics
through `sanntis.train.metrics` — so every number here is computed by the same code users run.

Each objective is a directory of numbered bash runbooks (`NN_*.sh`) that call tools directly, plus a few small Python
workers where no tool does the job. Every input is downloaded from a public repository by accession (MIBiG, NCBI,
GTDB) and checked against a sha256 or a recorded count, and every result is computed from those inputs on an empty
`DATA_ROOT`, InterProScan included. Expected results sit next to the runbooks.

## Objectives

| | Objective | Produces |
|---|---|---|
| [00](objectives/00-data-collection/) | Data collection | MIBiG archives + loci, taxonomy, negative-pool genomes, test backgrounds |
| [01](objectives/01-production-models/) | Production models | negative pool v2 (lineage-matched) and the SanntiS 1.2.0 model bundle (MIBiG 5.0rc1 active + pending as MIBiG-format archives) with its validation report |
| [02](objectives/02-mibig31-validation/) | MIBiG 3.1 validation | a model trained on v3-active, scored on Test A / B with isolate and gene-shuffled backgrounds (incl. building those sets) |
| [03](objectives/03-tool-comparison/) | Tool comparison | detection (antiSMASH / fungiSMASH / plantiSMASH, GECCO, SanntiS 0.9.5, BGC-Prophet) and classification (BGCat) against SanntiS 1.2.0: tables and figures |

Each objective's README lists its inputs (with accessions), steps, outputs, expected numbers and resources.

## Data and models

| | Where | Licence |
|---|---|---|
| SanntiS software | [github.com/Finn-Lab/SanntiS](https://github.com/Finn-Lab/SanntiS), release `v1.2.0`; Zenodo [10.5281/zenodo.24003727](https://doi.org/10.5281/zenodo.24003727) | Apache-2.0 |
| SanntiS 1.2.0 model | [huggingface.co/microbiome-informatics/sanntis](https://huggingface.co/microbiome-informatics/sanntis), tag `v1.2.0`; Zenodo `sanntis_model` | Apache-2.0 |
| Training data | [huggingface.co/datasets/microbiome-informatics/sanntis-training-data](https://huggingface.co/datasets/microbiome-informatics/sanntis-training-data), tags `negatives-v2`, `npclassifier-v1`, `mibig-5.0rc1-train`; Zenodo `mibig5_train` | CC BY 4.0 (`negatives/`: CC BY-SA 4.0) |
| Test sets, MIBiG 3.1 evaluation | Zenodo `nucleotide_testsets_v3`, `v3active_eval_v3` | CC BY 4.0 |

The Zenodo records (DOIs in the `record` column of [`env/zenodo.tsv`](env/zenodo.tsv)) are in the
[EUREMAP community](https://zenodo.org/communities/euremap). They hold what the runbooks produce; the runbooks do not
read them.

## Setup

```bash
bash env/install.sh                 # conda env + pinned SanntiS (env/sanntis.ref)
conda activate sanntisv1-analysis
export DATA_ROOT=/path/with/space   # downloads and results (default: <project area on NFS>/data, else ./data)
export RUNNER=local                 # or slurm: each step becomes `sbatch --wait`
export CONTAINER_ENGINE=docker      # or singularity / apptainer
export NF_CONFIG=...                # Nextflow config(s) for SanntiS's fan-out and InterProScan 6 (Codon: env/nextflow/codon.config)
bash objectives/00-data-collection/01_mibig.sh
```

Tool containers are pinned by digest in [`env/containers.tsv`](env/containers.tsv). Steps are idempotent: a finished
step leaves `$WORK/.done/<step>` and is skipped on re-run (`FORCE=1` re-runs it); logs go to `$WORK/logs/`.
With `RUNNER=slurm` each step is a job that re-runs the runbook for that step only, so `DATA_ROOT` must be on a
filesystem the compute nodes see; extra `sbatch` options go in `SLURM_ARGS`. `LIMIT=n` runs every objective as a
trial on a separate, empty `DATA_ROOT`: the first `n` genomes per kingdom, `10 x n` training BGCs, one detector
epoch, the first `n` test BGCs and contigs.

InterProScan 6 runs through SanntiS with 4 cores per task; with `RUNNER=slurm`, 16 chunk runs at a time with 40
queued tasks each (about 2,500 cores; `SANNTIS_IPS6_*` in `common/env.sh`). Its precalculated-match lookup is used
only when the InterPro Matches API serves the pinned release (InterPro 107.0); otherwise every protein is computed
locally.

## Layout

```
common/        env.sh (settings), lib.sh (fetch + checksum, NCBI / Zenodo download, containers, step runner)
env/           conda environment, SanntiS pin, container pins, Nextflow config (Codon)
objectives/    NN-<objective>/{README.md, NN_*.sh, workers/, expected/}
docs/          zenodo/ (deposit staging, metadata, NOTICE files)
```

## Data locations

On the EBI cluster everything lives in the project area `/nfs/production/rdf/metagenomics/projects/sanntisv1-data-and-analysis/`
(`PROJECT_ROOT`; `common/env.sh` picks it up when it exists and is writable, otherwise `DATA_ROOT` falls back to `./data`):

| Directory | Content |
|---|---|
| `data/` (`DATA_ROOT`) | everything the runbooks download or compute from public sources: `raw/`, `zenodo/`, `03/`, `cache/`, `work/`; can be deleted and rebuilt |

## Acknowledgements

This work is part of EUREMAP (European Research Infrastructure for Marine Bioprospecting, https://euremap.eu),
funded by the European Union's Horizon Europe research and innovation programme under grant agreement No 101131663
([doi:10.3030/101131663](https://doi.org/10.3030/101131663)). Views and opinions expressed are however those of the
author(s) only and do not necessarily reflect those of the European Union. Neither the European Union nor the granting
authority can be held responsible for them.

The analyses rely on public resources and their teams: MIBiG, NCBI / INSDC, GTDB, InterPro and NPClassifier (sources, licences and citations in `NOTICE`), and on the tools benchmarked in objective 03
(`env/containers.tsv`, [objective 03](objectives/03-tool-comparison/)). Computations ran on the EMBL-EBI Codon cluster.

## Licence

Code and the tables this repository writes: Apache-2.0 (`LICENSE`). Third-party content keeps its source's terms,
e.g. the GTDB r226 lineages in `testset_background_genomes.tsv` are CC BY-SA 4.0 (`NOTICE`). No third-party code,
image or model weights are redistributed; tools are referenced by version / tag, URL, sha256 and citation, and
fetched from their publishers at run time.
