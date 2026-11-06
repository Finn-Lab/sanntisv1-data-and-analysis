#!/usr/bin/env bash
# Genomes of negative pool v1 (205 GTDB prokaryotes, 23 RefSeq eukaryotes) from NCBI Datasets (md5-checked package).
#   → $DATA_ROOT/raw/ncbi/assemblies/<acc>/<acc>.fna.gz (+ .gff.gz, assembly_data_report.jsonl)
# Optional: objectives 01/02 download pool v1 itself from Zenodo. These genomes are what `sanntis build-negatives`
# rebuilds it from (pool v1 used the MGnify pipeline gene calls, so a rebuild is not byte-identical).
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need curl unzip md5sum
list=$REPO/objectives/00-data-collection/accessions/pool_v1_genomes.tsv
out=$DATA_ROOT/raw/ncbi/assemblies

genomes() {
    rows "$list" | awk -v o="$out" '{print $1, o}' | each "${JOBS:-4}" fetch_assembly
    local n
    n=$(rows "$list" | cut -f1 | while read -r a; do [[ -e $out/$a/done ]] && echo; done | wc -l)
    expect_count "pool v1 genomes" "$n" "$(rows "$list" | wc -l)"
}

step pool_genomes 4 4G 12:00:00 genomes
