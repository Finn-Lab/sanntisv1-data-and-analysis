# shellcheck shell=bash
# Source from a runbook:  source "$(dirname "$0")/../../common/env.sh"
# Every setting can be overridden from the environment before sourcing.

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export REPO
# shellcheck disable=SC2034  # used by step() in lib.sh
SCRIPT=$(realpath "$0")   # the runbook and its arguments (SLURM steps re-run it)
# shellcheck disable=SC2034
SCRIPT_ARGS=("$@")

PROJECT_ROOT=${PROJECT_ROOT:-/nfs/production/rdf/metagenomics/projects/sanntisv1-data-and-analysis}  # shared project area on NFS
# downloads + results (large; outside git): the project area when it exists, else ./data
if [[ -d $PROJECT_ROOT/data && -w $PROJECT_ROOT/data ]]; then DATA_ROOT=${DATA_ROOT:-$PROJECT_ROOT/data}; else DATA_ROOT=${DATA_ROOT:-$REPO/data}; fi
export DATA_ROOT
export WORK=${WORK:-$DATA_ROOT/work}                 # step markers + logs
export RUNNER=${RUNNER:-local}                       # local | slurm
export CONTAINER_ENGINE=${CONTAINER_ENGINE:-docker}  # docker | singularity | apptainer
export SANNTIS_CACHE=${SANNTIS_CACHE:-$DATA_ROOT/cache/sanntis}
export PYTHON=${PYTHON:-python}                      # the analysis environment's interpreter (env/install.sh)
export PYTHONPATH=$REPO/workers${PYTHONPATH:+:$PYTHONPATH}
export NF_CONFIG=${NF_CONFIG:-}                      # Nextflow config(s) for SanntiS / InterProScan 6 (space-separated)

set -euo pipefail
# shellcheck source=SCRIPTDIR/lib.sh
source "$REPO/common/lib.sh"
mkdir -p "$DATA_ROOT" "$WORK"
