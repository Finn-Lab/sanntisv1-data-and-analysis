#!/usr/bin/env bash
# Inputs: objective 02's nucleotide test sets (deposit nucleotide_testsets_v2) and its v3-active evaluation (the
# model bundle, test-BGC table, NPC labels, SanntiS tokens, scores): $DATA_ROOT/02/eval when this DATA_ROOT holds
# an objective 02 run, else V3ACTIVE_EVAL=<such a directory>, else the v3active_eval_v2 deposit.
# Tool outputs are not fetched: 02_tools.sh runs the tools.
#   → $DATA_ROOT/03/inputs/{nucleotide_testsets,v3active_eval}/
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need curl tar sha256sum
in=$DATA_ROOT/03/inputs

inputs() {
    deposit nucleotide_testsets_v2 "$in"
    local ev=${V3ACTIVE_EVAL:-$DATA_ROOT/02/eval}
    if [[ -d $ev/bundle && -s $ev/test_bgcs.parquet ]]; then
        log "v3active_eval: $ev"
        ln -sfn "$(realpath "$ev")" "$in/v3active_eval"
    else
        [[ -z ${V3ACTIVE_EVAL:-} ]] || die "V3ACTIVE_EVAL: no bundle/ + test_bgcs.parquet under $V3ACTIVE_EVAL"
        deposit v3active_eval_v2 "$in"
    fi
    expect_count "test contigs" "$(find "$in/nucleotide_testsets" -name '*.gbk' | wc -l)" 445
    expect_count "test BGCs" "$("$PYTHON" -c 'import pandas as pd, sys; print(len(pd.read_parquet(sys.argv[1])))' \
        "$in/v3active_eval/test_bgcs.parquet")" 375
}

step 03_inputs 1 4G 2:00:00 inputs
