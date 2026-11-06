#!/usr/bin/env bash
# Training positives: MIBiG 5.0rc1 active + pending as a MIBiG-format archive pair for `sanntis train --external-gbk /
# --external-json` (5.0rc1 has no GBK archive, so `--mibig 5.0rc1` cannot be used). Built from objective 00's MIBiG
# downloads (default), or fetched frozen from Zenodo (MIBIG_TRAIN=deposit). Either way checked against expected/.
#   → $DATA_ROOT/01/inputs/mibig_5.0rc1-train/{mibig_5.0rc1-train_{gbk,json}.tar.gz, manifest.tsv, excluded.tsv,
#     provenance.json}
#   → $DATA_ROOT/01/deposit/sanntisv1_mibig5.0rc1_train.tar.gz (a build only: the Zenodo file)
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need tar gzip sha256sum
here=$REPO/objectives/01-production-models
acc=$REPO/objectives/00-data-collection/accessions
in=$DATA_ROOT/01/inputs
bgcs=$in/mibig_5.0rc1-train  # not `name`: step() has a local of that name

check() {  # the archives hold what expected/ lists, BGC by BGC (source, locus, length, sequence and protein md5)
    diff -q "$here/expected/mibig_train_manifest.tsv" "$bgcs/manifest.tsv" >/dev/null \
        || die "$bgcs/manifest.tsv differs from expected/mibig_train_manifest.tsv"
    diff -q "$here/expected/mibig_train_excluded.tsv" "$bgcs/excluded.tsv" >/dev/null \
        || die "$bgcs/excluded.tsv differs from expected/mibig_train_excluded.tsv"
    expect_count "MIBiG 5.0rc1 training BGCs" "$(rows "$bgcs/manifest.tsv" | wc -l)" 2795
}

build() {
    local raw=$DATA_ROOT/raw/mibig dep=$DATA_ROOT/01/deposit/sanntisv1_mibig5.0rc1_train.tar.gz sha want
    [[ -s $raw/mibig_json_5.0rc1.tar.bz2 && -s $raw/mibig_gbk_4.0.tar.gz ]] \
        || die "no MIBiG archives in $raw (run objectives/00-data-collection/01_mibig.sh)"
    expect_count "objective-00 nuccore loci" "$(find "$raw/nuccore" -name '*.gbk' | wc -l)" \
        "$(awk 'NR > 1' "$acc/mibig_nuccore_loci.tsv" | cut -f2- | sort -u | wc -l)"
    rm -rf "$bgcs"
    "$PYTHON" "$here/workers/mibig_train_archives.py" "$raw" "$bgcs"
    check
    # the Zenodo file: deterministic, so a build with the same SanntiS / Biopython reproduces its sha256
    mkdir -p "$(dirname "$dep")"
    tar -C "$in" --sort=name --mtime=@0 --owner=0 --group=0 --numeric-owner -cf - "$(basename "$bgcs")" | gzip -n > "$dep"
    sha=$(sha256sum "$dep" | cut -d' ' -f1)
    want=$(awk -F'\t' '$1 == "mibig5_train" {print $4}' "$REPO/env/zenodo.tsv")
    if [[ $sha == "$want" ]]; then log "deposit sha256 = env/zenodo.tsv ✓"
    else log "deposit sha256 $sha != env/zenodo.tsv ($want): other SanntiS / Biopython? content matches expected/"; fi
}

fetch() {
    deposit mibig5_train "$in"
    check
}

case ${MIBIG_TRAIN:-build} in
    build)   step 01_mibig_train 2 8G 1:00:00 build ;;
    deposit) step 01_mibig_train 1 2G 1:00:00 fetch ;;
    *) die "MIBIG_TRAIN must be build or deposit (got '$MIBIG_TRAIN')" ;;
esac
