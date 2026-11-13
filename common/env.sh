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
# LIMIT=n (a trial on a separate DATA_ROOT): first rows of the tables, 10 x n training BGCs, 1 detector epoch
[[ -z ${LIMIT:-} ]] || export SANNTIS_DEBUG_MAX_EPOCHS=${SANNTIS_DEBUG_MAX_EPOCHS:-1}
# InterProScan 6 through SanntiS: 4 cores per task; on SLURM 16 chunk runs at once x 40 queued tasks each (about
# 2,500 cores; one Nextflow head per chunk run in the step's job, see the step memory); the precalculated-match
# lookup only when the Matches API serves the pinned InterPro release (else every protein is computed locally)
export SANNTIS_IPS6_CPUS=${SANNTIS_IPS6_CPUS:-4}
export SANNTIS_IPS6_LOOKUP=${SANNTIS_IPS6_LOOKUP:-auto}
if [[ $RUNNER == slurm ]]; then
    export SANNTIS_IPS6_PARALLEL=${SANNTIS_IPS6_PARALLEL:-16} SANNTIS_IPS6_QUEUE_SIZE=${SANNTIS_IPS6_QUEUE_SIZE:-40}
    export SANNTIS_IPS6_SUBMIT_RATE=${SANNTIS_IPS6_SUBMIT_RATE:-100/1min}
    export NXF_OPTS=${NXF_OPTS:--Xms256m -Xmx2g}
fi

set -euo pipefail
# shellcheck source=SCRIPTDIR/lib.sh
source "$REPO/common/lib.sh"
mkdir -p "$DATA_ROOT" "$WORK"
