#!/usr/bin/env bash
# The production bundle: `sanntis train` on 01_inputs.sh's MIBiG 5.0rc1 archives (external data, `--mibig none`) and
# negative pool v2 (03_pool.sh), default protocol (5-split CV → epochs and NPC thresholds, t* = 0.20 from the recipe →
# 3 seeds → bundle + out-of-fold report). InterProScan 6 on the positives is cached in SANNTIS_CACHE by protein md5
# (02_plan.sh computed it). Fan-out units (CV splits, seeds; label `train`) and InterProScan 6 go through Nextflow
# with NF_CONFIG; with RUNNER=slurm this step is the orchestrator job (small, long).
#   NEGATIVES=<pool dir | Hub revision>   (default $DATA_ROOT/01/pool/v2)
#   TRAIN_ARGS="..."                      extra `sanntis train` options, e.g. "--runner inline --workers 4", "--gpu"
#   → $DATA_ROOT/01/run/{bundle/, report/report.md, tune/tune.json}
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need sanntis nextflow
here=$REPO/objectives/01-production-models
d=$DATA_ROOT/01/inputs/mibig_5.0rc1-train
run=$DATA_ROOT/01/run

train() {
    [[ -s $d/mibig_5.0rc1-train_gbk.tar.gz ]] || die "no training archives in $d (run 01_inputs.sh)"
    mapfile -t nf < <(nf_args)
    mapfile -t pos < <(train_archives "$d" mibig_5.0rc1-train)
    # shellcheck disable=SC2086  # TRAIN_ARGS: options, split on purpose
    SANNTIS_THREADS=${SANNTIS_THREADS:-4} sanntis -v train --run "$run" --mibig none \
        --external-gbk "${pos[0]}" --external-json "${pos[1]}" \
        --negatives "${NEGATIVES:-$DATA_ROOT/01/pool/v2}" --container-engine "$CONTAINER_ENGINE" "${nf[@]}" ${TRAIN_ARGS:-}
    [[ -s $run/report/report.md ]] || die "no report in $run/report"
    [[ -z ${LIMIT:-} ]] || return 0
    # collect: 2,795 archived BGCs minus the 18 loci > 200 kb
    expect_count "collected loci" "$("$PYTHON" -c 'import json, sys; print(json.load(open(sys.argv[1]))["counts"]["loci"])' \
        "$run/collect/provenance.json")" 2777
    compare_expected "$run/bundle/report" "$here/expected" validation.tsv validation_strata.tsv \
        validation_class_curves.tsv sweep.tsv npc_validation.tsv
}

step 01_train 4 48G 7-00:00:00 train
