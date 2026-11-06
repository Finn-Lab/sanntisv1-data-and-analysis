#!/usr/bin/env bash
# Create the analysis environment and install the pinned SanntiS:  bash env/install.sh [env name]
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
name=${1:-sanntisv1-analysis}
read -r ref url < <(grep -v '^#' "$here/sanntis.ref")
conda env create -n "$name" -f "$here/environment.yml"
conda run -n "$name" pip install "sanntis[train] @ git+$url@$ref"
echo "conda activate $name"
