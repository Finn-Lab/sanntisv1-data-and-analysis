#!/usr/bin/env bash
# MIBiG archives (sha256-pinned) + the loci missing from the 4.0 GBK archive (NCBI nuccore).
#   → $DATA_ROOT/raw/mibig/<archive>, $DATA_ROOT/raw/mibig/nuccore/<acc>_<from>_<to>.gbk
# SanntiS fetches the same files itself (`sanntis train --mibig <release>`); objectives 01/02 point SANNTIS_CACHE here.
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need curl sha256sum
acc=$REPO/objectives/00-data-collection/accessions
out=$DATA_ROOT/raw/mibig

archives() {
    rows "$acc/mibig_archives.tsv" | while IFS=$'\t' read -r _ _ url sha; do
        fetch "$url" "$out/$(basename "$url")" "$sha"
    done
}

nuccore() {
    rows "$acc/mibig_nuccore_loci.tsv" | while IFS=$'\t' read -r _ nuc from to; do
        fetch_nuccore "$nuc" "$from" "$to" "$out/nuccore/${nuc}_${from}_${to}.gbk"
    done
    # distinct loci: two entries can share one (BGC0003257 / BGC0003262)
    expect_count "nuccore loci" "$(find "$out/nuccore" -name '*.gbk' | wc -l)" \
        "$(rows "$acc/mibig_nuccore_loci.tsv" | cut -f2- | sort -u | wc -l)"
    sums "$out/nuccore" '*.gbk'
}

step mibig_archives 1 2G 1:00:00 archives
step mibig_nuccore 1 2G 4:00:00 nuccore
