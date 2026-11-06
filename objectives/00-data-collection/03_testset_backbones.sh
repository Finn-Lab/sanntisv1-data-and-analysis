#!/usr/bin/env bash
# Backbone genomes of the nucleotide test sets (v5 isolate: 55 proka + 15 euk; v5shuf2: 102 + 22; 138 distinct).
#   → $DATA_ROOT/raw/ncbi/assemblies/<acc>/ (shared with 02_pool_genomes.sh)
# Optional: objective 03 downloads the built test sets from Zenodo. These genomes are what the builders start from.
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need curl unzip md5sum
list=$REPO/objectives/00-data-collection/accessions/testsets_nucleotide_backbones.tsv
out=$DATA_ROOT/raw/ncbi/assemblies

backbones() {
    local accs n
    accs=$(rows "$list" | cut -f3 | sort -u)
    awk -v o="$out" '{print $1, o}' <<< "$accs" | each "${JOBS:-4}" fetch_assembly
    n=$(while read -r a; do [[ -e $out/$a/done ]] && echo; done <<< "$accs" | wc -l)
    expect_count "test-set backbone genomes" "$n" "$(wc -l <<< "$accs")"
}

step testset_backbones 4 4G 12:00:00 backbones
