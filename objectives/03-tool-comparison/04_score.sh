#!/usr/bin/env bash
# Score SanntiS 1.2.0 (v3-active model at its t*) and every tool run on the same contigs and CDS, then compare with
# expected/ once its tables exist for these sets (expected/v1sets/ holds the previous sets' tables).
#   → $DATA_ROOT/03/score/{results.pkl,summary.tsv,metrics.tsv,strata.tsv,radar.tsv,curve.tsv,auprc.tsv,paired.tsv}
#   → $DATA_ROOT/03/score/paired_bgc.tsv: SanntiS 1.2.0 − each tool, paired bootstrap of the BGC-level metrics on
#     B-shuffled (scopes all / proka / euk / fungi)
# Tool outputs: the runs of 02_tools.sh ($DATA_ROOT/03/tools); a tool without a run is left out (SanntiS 1.2.0 is
# always scored).
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
o=$REPO/objectives/03-tool-comparison
in=$DATA_ROOT/03/inputs
out=$DATA_ROOT/03/score

inputs() {  # the arguments shared by score.py and paired_bgc.py
    local d
    for d in "$DATA_ROOT/03/tools"; do [[ ! -d $d ]] || ARGS+=(--tools "$d"); done
    ARGS+=(--testsets "$in/nucleotide_testsets" --model "$in/v3active_eval/bundle"
           --test-bgcs "$in/v3active_eval/test_bgcs.parquet" --tokens "$DATA_ROOT/03/sanntis/tokens.parquet")
}

score() {
    local ARGS=() files=()
    inputs
    [[ " ${ARGS[*]} " == *" --tools "* ]] || log "no tool runs under $DATA_ROOT/03/tools (02_tools.sh): scoring SanntiS 1.2.0 only"
    "$PYTHON" "$o/workers/score.py" "${ARGS[@]}" --out "$out"
    if [[ -s $o/expected/summary.tsv ]]; then
        # paired_bgc.tsv is compared by its own step
        mapfile -t files < <(cd "$o/expected" && ls -- *.tsv | grep -vx paired_bgc.tsv)
        "$PYTHON" -m sv1.expected compare "$out" "$o/expected" --files "${files[@]}"
    else
        log "expected/ not frozen for these sets yet: no comparison"
    fi
}

paired_bgc() {
    local ARGS=()
    inputs
    "$PYTHON" "$o/workers/paired_bgc.py" "${ARGS[@]}" --out "$out/paired_bgc.tsv"
    if [[ -s $o/expected/paired_bgc.tsv ]]; then
        "$PYTHON" -m sv1.expected compare "$out" "$o/expected" --files paired_bgc.tsv
    else
        log "expected/paired_bgc.tsv not frozen yet: no comparison"
    fi
}

step 03_score 4 16G 6:00:00 score
step 03_paired_bgc 4 16G 1:00:00 paired_bgc
