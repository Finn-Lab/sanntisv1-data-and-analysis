#!/usr/bin/env bash
# SanntiS tokens of the test-set proteins: objective 02's (InterProScan 6.0.3.0 / InterPro 107.0 through SanntiS,
# 05_score.sh; md5-keyed).
#   → $DATA_ROOT/03/sanntis/tokens.parquet
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
in=$DATA_ROOT/03/inputs
out=$DATA_ROOT/03/sanntis

tokens() {
    mkdir -p "$out"
    [[ -s $in/v3active_eval/tokens.parquet ]] || die "no tokens in $in/v3active_eval (objective 02's 05_score.sh)"
    cp -f "$in/v3active_eval/tokens.parquet" "$out/tokens.parquet"
}

step 03_sanntis_tokens 2 16G 24:00:00 tokens
