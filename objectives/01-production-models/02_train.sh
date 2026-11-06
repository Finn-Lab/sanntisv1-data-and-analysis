#!/usr/bin/env bash
# The production bundle: `sanntis train` on 01_inputs.sh's MIBiG 5.0rc1 archives (external data, `--mibig none`) and
# negative pool v1, default protocol (5-split CV tunes epochs / t* / NPC thresholds → 3 seeds → bundle + out-of-fold
# report). Needs `sanntis setup` (InterProScan 6 data) once; InterProScan 6 on the positives is cached in
# SANNTIS_CACHE by protein md5. Fan-out units (CV splits, seeds; label `train`) and InterProScan 6 go through
# Nextflow with NF_CONFIG; with RUNNER=slurm this step is the orchestrator job (small, long).
#   NEGATIVES=<pool dir | Hub revision>   (default v1)
#   TRAIN_ARGS="..."                      extra `sanntis train` options, e.g. "--runner inline --workers 4", "--gpu"
#   → $DATA_ROOT/01/run/{bundle/, report/report.md, tune/tune.json}
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need sanntis nextflow
d=$DATA_ROOT/01/inputs/mibig_5.0rc1-train
run=$DATA_ROOT/01/run

train() {
    [[ -s $d/mibig_5.0rc1-train_gbk.tar.gz ]] || die "no training archives in $d (run 01_inputs.sh)"
    mapfile -t nf < <(nf_args)
    # shellcheck disable=SC2086  # TRAIN_ARGS: options, split on purpose
    SANNTIS_THREADS=${SANNTIS_THREADS:-4} sanntis -v train --run "$run" --mibig none \
        --external-gbk "$d/mibig_5.0rc1-train_gbk.tar.gz" --external-json "$d/mibig_5.0rc1-train_json.tar.gz" \
        --negatives "${NEGATIVES:-v1}" --container-engine "$CONTAINER_ENGINE" "${nf[@]}" ${TRAIN_ARGS:-}
    # collect = the release pool: 2,795 archived BGCs minus the 18 loci > 200 kb
    expect_count "collected loci" "$("$PYTHON" -c 'import json, sys; print(json.load(open(sys.argv[1]))["counts"]["loci"])' \
        "$run/collect/provenance.json")" 2777
    [[ -s $run/report/report.md ]] || die "no report in $run/report"
}

step 01_train 4 48G 7-00:00:00 train
