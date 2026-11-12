#!/usr/bin/env bash
# Product classification of the test BGCs: SanntiS 1.0 NPC (v3-active model, true loci) vs BGCat (02_tools.sh bgcat
# on the v2 isolate set), in the SanntiS NPC label space.
#   → $DATA_ROOT/03/classification/{classification.tsv,classification_per_bgc.tsv}
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
o=$REPO/objectives/03-tool-comparison
in=$DATA_ROOT/03/inputs
out=$DATA_ROOT/03/classification
# NPClassifier ontology (class → superclass → pathway), MIT, pinned to a commit
NPC_INDEX_URL=https://raw.githubusercontent.com/mwang87/NP-Classifier/78e52f1d27484841b303417dc2847430036092d6/Classifier/dict/index_v1.json
NPC_INDEX_SHA=4ed2a982b6ef75372c5c4e69722f3f4ab239c8495d6cff2a6dec4fa95605d7cd

classification() {
    local pred=$DATA_ROOT/03/tools/bgcat_4f77517/predictions.json
    [[ -s $pred ]] || die "no BGCat predictions in $pred (run 02_tools.sh antismash-relaxed, then bgcat)"
    fetch "$NPC_INDEX_URL" "$DATA_ROOT/cache/npclassifier/index_v1.json" "$NPC_INDEX_SHA"
    # MIBiG 4.0 accessions: BGCat's released weights may have seen them, so a subset without them is scored too
    local v4 url sha
    read -r url sha < <(awk -F'\t' '$1 == "4.0" && $2 == "json" {print $3, $4}' \
        "$REPO/objectives/00-data-collection/accessions/mibig_archives.tsv")
    v4=$DATA_ROOT/raw/mibig/$(basename "$url")
    fetch "$url" "$v4" "$sha"
    "$PYTHON" "$o/workers/npc_compare.py" --test-bgcs "$in/v3active_eval/test_bgcs.parquet" \
        --labels "$in/v3active_eval/labels" --model "$in/v3active_eval/bundle" --bgcat "$pred" \
        --npc-index "$DATA_ROOT/cache/npclassifier/index_v1.json" --mibig4-json "$v4" --out "$out"
    # the frozen tables (expected/classification/, from these BGCat predictions; PROVENANCE.json)
    [[ -s $out/classification.tsv ]] || { [[ -n ${LIMIT:-} ]] && return 0; die "no $out/classification.tsv"; }
    compare_expected "$out" "$o/expected/classification" classification.tsv
}

step 03_classification 2 8G 1:00:00 classification
