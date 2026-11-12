#!/usr/bin/env bash
# Genomes of negative pool v2 (accessions/pool_v2_genomes.tsv, written by `sanntis plan-negatives`; objective 01
# re-derives it): prokaryotes as NCBI genome FASTA (genes are called by `sanntis build-negatives`), fungi with their
# annotated RefSeq GenBank; NCBI Datasets, md5-checked package.
#   → $DATA_ROOT/raw/ncbi/assemblies/<acc>/<acc>.fna.gz (+ .gff.gz, euk: .gbff.gz, assembly_data_report.jsonl)
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need curl unzip md5sum
list=$REPO/objectives/00-data-collection/accessions/pool_v2_genomes.tsv
out=$DATA_ROOT/raw/ncbi/assemblies

genomes() {
    local n
    # columns: genome kingdom lineage source accession ...; no trailing blank field (`xargs -L 1` would join lines)
    rows_per_kingdom "$list" 2 | awk -F'\t' -v o="$out" '$2 == "euk" {print $5, o, "gbff"; next} {print $5, o}' \
        | each "${JOBS:-4}" fetch_assembly
    n=$(rows_per_kingdom "$list" 2 | cut -f5 | while read -r a; do [[ -e $out/$a/done ]] && echo; done | wc -l)
    expect_count "pool v2 genomes" "$n" "$(rows_per_kingdom "$list" 2 | wc -l)"
}

step pool_genomes 4 4G 24:00:00 genomes
