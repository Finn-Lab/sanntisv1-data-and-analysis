#!/usr/bin/env bash
# The v3-active model: `sanntis train` on 01_v3active.sh's archives (external data, `--mibig none`) and negative pool
# v1, default protocol (5-split CV tunes epochs / t* / NPC thresholds → 3 seeds → bundle + out-of-fold report).
# Needs `sanntis setup` (InterProScan 6 data) once. Runs as objective 01's 02_train.sh (NF_CONFIG, RUNNER=slurm).
#   NEGATIVES=<pool dir | Hub revision>   (default v1)
#   TRAIN_ARGS="..."                      extra `sanntis train` options
#   → $DATA_ROOT/02/run/{bundle/, report/report.md, tune/tune.json}
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need sanntis nextflow
v3=$DATA_ROOT/02/inputs/mibig_v3-active
run=$DATA_ROOT/02/run

train() {
    [[ -s $v3/mibig_v3-active_gbk.tar.gz ]] || die "no v3-active archives in $v3 (run 01_v3active.sh)"
    mapfile -t nf < <(nf_args)
    # shellcheck disable=SC2086  # TRAIN_ARGS: options, split on purpose
    SANNTIS_THREADS=${SANNTIS_THREADS:-4} sanntis -v train --run "$run" --mibig none \
        --external-gbk "$v3/mibig_v3-active_gbk.tar.gz" --external-json "$v3/mibig_v3-active_json.tar.gz" \
        --negatives "${NEGATIVES:-v1}" --container-engine "$CONTAINER_ENGINE" "${nf[@]}" ${TRAIN_ARGS:-}
    # 2,419 archived BGCs minus the 17 loci > 200 kb
    expect_count "collected loci" "$("$PYTHON" -c 'import json, sys; print(json.load(open(sys.argv[1]))["counts"]["loci"])' \
        "$run/collect/provenance.json")" 2402
    [[ -s $run/report/report.md ]] || die "no report in $run/report"
}

step 02_train 4 48G 7-00:00:00 train
