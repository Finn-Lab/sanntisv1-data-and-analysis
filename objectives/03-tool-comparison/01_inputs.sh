#!/usr/bin/env bash
# Inputs: objective 02's nucleotide test sets and its v3-active evaluation (model bundle, test-BGC table, NPC labels,
# SanntiS tokens, scores) from this DATA_ROOT. Tool outputs are not fetched: 02_tools.sh runs the tools.
#   → $DATA_ROOT/03/inputs/{nucleotide_testsets,v3active_eval} (links)
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
in=$DATA_ROOT/03/inputs

inputs() {
    local ts=$DATA_ROOT/02/testsets/testsets ev=$DATA_ROOT/02/eval
    [[ -s $ts/README.txt ]] || die "no test sets in $ts (run objectives/02-mibig31-validation/04_testsets.sh)"
    [[ -d $ev/bundle && -s $ev/test_bgcs.parquet ]] || die "no evaluation in $ev (run objectives/02-mibig31-validation/05_score.sh)"
    mkdir -p "$in"
    ln -sfn "$(realpath "$ts")" "$in/nucleotide_testsets"
    ln -sfn "$(realpath "$ev")" "$in/v3active_eval"
    [[ -z ${LIMIT:-} ]] || return 0
    expect_count "test contigs" "$(find "$in/nucleotide_testsets/" -name '*.gbk' | wc -l)" 445
    expect_count "test BGCs" "$("$PYTHON" -c 'import pandas as pd, sys; print(len(pd.read_parquet(sys.argv[1])))' \
        "$in/v3active_eval/test_bgcs.parquet")" 375
}

step 03_inputs 1 4G 2:00:00 inputs
