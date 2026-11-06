#!/usr/bin/env bash
# SanntiS tokens of the test-set proteins: objective 02 already annotated them (InterProScan 6.0.3.0 / InterPro
# 107.0 through SanntiS, 05_score.sh) and ships them in v3active_eval/tokens.parquet (md5-keyed, so they apply to
# any copy of the contigs); without that file they are computed here (~96k proteins, cached in SANNTIS_CACHE).
#   → $DATA_ROOT/03/sanntis/tokens.parquet
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
w=$REPO/objectives/03-tool-comparison/workers
in=$DATA_ROOT/03/inputs
out=$DATA_ROOT/03/sanntis

tokens() {
    mkdir -p "$out"
    if [[ -s $in/v3active_eval/tokens.parquet ]]; then
        log "tokens from v3active_eval"
        cp -f "$in/v3active_eval/tokens.parquet" "$out/tokens.parquet"
        return 0
    fi
    need nextflow
    mapfile -t nf < <(nf_args)
    "$PYTHON" "$w/annotate_ips6.py" --testsets "$in/nucleotide_testsets" --out "$out/tokens.parquet" \
        --container-engine "$CONTAINER_ENGINE" "${nf[@]}"
}

step 03_sanntis_tokens 2 16G 24:00:00 tokens
