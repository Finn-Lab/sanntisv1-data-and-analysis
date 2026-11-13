#!/usr/bin/env bash
# Score the v3-active model on the nucleotide test sets (04_testsets.sh): first the SanntiS tokens of every test
# protein (InterProScan 6 through SanntiS; fan-out via Nextflow with NF_CONFIG, cached in SANNTIS_CACHE), then
# Test A / B on the isolate and shuffled contigs at the bundle's t* (`sanntis.train.metrics`: metrics with 95 % CIs
# over source contigs, strata, radar, PR curves, NPC classification), and the evaluation bundle objective 03 reads.
#   MODEL=<bundle dir>   (default $DATA_ROOT/02/run/bundle, the 02_train.sh model)
#   TRAIN_ENTITIES=<entities.parquet>  the training positives whose Pfam sets define novelty (default: the 02 run)
#   → $DATA_ROOT/02/eval/{tokens.parquet, ips6_run/, bundle -> MODEL, test_bgcs.parquet, labels/, score/}
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need nextflow
here=$REPO/objectives/02-mibig31-validation
ts=$DATA_ROOT/02/testsets
eval=$DATA_ROOT/02/eval

tokens() {
    [[ -s $ts/test_bgcs.parquet ]] || die "no test sets in $ts (run 04_testsets.sh)"
    mapfile -t nf < <(nf_args)
    "$PYTHON" "$here/workers/tokens.py" --testsets "$ts/testsets" --out "$eval/tokens.parquet" \
        --container-engine "$CONTAINER_ENGINE" "${nf[@]}"
}

step 02_tokens 2 48G 2-00:00:00 tokens  # memory: the InterProScan 6 Nextflow heads (common/env.sh)

score() {
    local model=${MODEL:-$DATA_ROOT/02/run/bundle}
    [[ -s $model/config.json ]] || die "no model bundle in $model (02_train.sh not finished?); MODEL=<dir> scores another"
    [[ -s $eval/tokens.parquet ]] || die "no tokens in $eval (02_tokens step)"
    [[ -s $DATA_ROOT/01/run/collect/loci.parquet ]] \
        || die "no objective-01 collect in $DATA_ROOT/01/run (the test BGCs' NPC labels; run its 02_train.sh to that step)"
    "$PYTHON" "$here/workers/score_testsets.py" --testsets "$ts/testsets" --test-bgcs "$ts/test_bgcs.parquet" \
        --tokens "$eval/tokens.parquet" --model "$model" \
        --train-entities "${TRAIN_ENTITIES:-$DATA_ROOT/02/run/entities/entities.parquet}" \
        --labels-loci "$DATA_ROOT/01/run/collect/loci.parquet" --out "$eval"
    ln -sfn "$(realpath "$model")" "$eval/bundle"  # objective 03: V3ACTIVE_EVAL=$eval
    expect_count "scored test sets" "$(find "$eval/score" -name '*.pkl' | wc -l)" 4
    column -t "$eval/score/summary.tsv"
    compare_expected "$eval/score" "$here/expected" test_auprc.tsv test_curve.tsv test_metrics.tsv \
        test_npc_per_pathway.tsv test_npc.tsv test_radar.tsv test_strata.tsv
}

step 02_score 4 16G 12:00:00 score
