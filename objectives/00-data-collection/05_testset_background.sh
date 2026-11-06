#!/usr/bin/env bash
# Background genomes of objective 02's test sets, none in negative pool v1 (102 GTDB r226 prokaryotes, 22 RefSeq
# eukaryotes; selection: objective 02 workers/sample_background.py) from NCBI Datasets (md5-checked package).
# Eukaryotes also with their annotated GenBank (RefSeq CDS); prokaryotes get gene calls in objective 02.
#   → $DATA_ROOT/raw/ncbi/assemblies/<acc>/<acc>.fna.gz (+ .gff.gz, euk: .gbff.gz, assembly_data_report.jsonl)
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need curl unzip md5sum
list=$REPO/objectives/00-data-collection/accessions/testset_background_genomes.tsv
out=$DATA_ROOT/raw/ncbi/assemblies

genomes() {
    local n
    # no trailing blank field: `xargs -L 1` would join the next line to it
    rows "$list" | awk -F'\t' -v o="$out" '$3 == "euk" {print $1, o, "gbff"; next} {print $1, o}' \
        | each "${JOBS:-4}" fetch_assembly
    n=$(rows "$list" | cut -f1 | while read -r a; do [[ -e $out/$a/done ]] && echo; done | wc -l)
    expect_count "test-background genomes" "$n" "$(rows "$list" | wc -l)"
}

step testset_background 4 4G 12:00:00 genomes
