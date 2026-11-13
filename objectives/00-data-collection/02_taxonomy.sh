#!/usr/bin/env bash
# Taxonomy of the negative-pool plan (objective 01) and of the test-background sampling: one NCBI taxdump of the
# NCBI taxdump archive and the GTDB r226 metadata tables, sha256-pinned in accessions/taxonomy.tsv.
#   → $DATA_ROOT/raw/ncbi/taxdump/taxdmp_<date>.zip, $SANNTIS_CACHE/gtdb/r226/*_metadata_r226.tsv.gz (where SanntiS
#     reads them)
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need curl sha256sum
LIMIT=  # three files, always complete
list=$REPO/objectives/00-data-collection/accessions/taxonomy.tsv

taxonomy() {
    local dest
    rows "$list" | while IFS=$'\t' read -r name url sha rel; do
        case $rel in
            cache/sanntis/*) dest=$SANNTIS_CACHE/${rel#cache/sanntis/} ;;
            *) dest=$DATA_ROOT/$rel ;;
        esac
        log "$name"
        fetch "$url" "$dest" "$sha"
    done
}

step taxonomy 1 2G 2:00:00 taxonomy
