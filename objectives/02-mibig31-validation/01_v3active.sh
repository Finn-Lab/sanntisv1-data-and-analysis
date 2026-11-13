#!/usr/bin/env bash
# v3-active training data: the MIBiG 3.1-era entries (accession < BGC0002750) of objective 01's 5.0rc1 archives, as a
# MIBiG-format archive pair for `sanntis train --external-gbk / --external-json`; checked against expected/.
#   → $DATA_ROOT/02/inputs/mibig_v3-active/{mibig_v3-active_{gbk,json}.tar.gz, manifest.tsv, provenance.json}
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
here=$REPO/objectives/02-mibig31-validation
src=$DATA_ROOT/01/inputs/mibig_5.0rc1-train
v3=$DATA_ROOT/02/inputs/mibig_v3-active

v3active() {
    [[ -s $src/manifest.tsv ]] || die "no objective-01 archives in $src (run objectives/01-production-models/01_inputs.sh)"
    rm -rf "$v3"
    "$PYTHON" "$here/workers/v3_archives.py" "$src" "$v3"
    diff -q "$here/expected/v3active_manifest.tsv" "$v3/manifest.tsv" >/dev/null \
        || die "$v3/manifest.tsv differs from expected/v3active_manifest.tsv"
    expect_count "v3-active BGCs" "$(awk 'NR > 1' "$v3/manifest.tsv" | wc -l)" 2419
}

step 02_v3active 1 4G 1:00:00 v3active
