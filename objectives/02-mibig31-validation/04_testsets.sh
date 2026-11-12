#!/usr/bin/env bash
# The nucleotide test sets: Test B (MIBiG 5.0rc1 active + pending >= BGC0002750, <= 200 kb: 375 BGCs; Test A = its
# 178 active entries) inserted into the 70 isolate test backbones of negative pool v2 (objective 01's 03_pool.sh;
# held out of training) as isolate contigs, and into gene-shuffled contigs of the test-background genomes
# (03_background.sh; none in the pool), one BGC each; worker build_testsets.py. Then the Zenodo deposit
# (deterministic tarball, sha256 pinned in env/zenodo.tsv as `nucleotide_testsets_v3`).
#   → $DATA_ROOT/02/testsets/{testsets/<set>/<track>/*.{gbk,cds.parquet,ground_truth.bed}, test_bgcs.parquet,
#     contigs.tsv, manifest.json}, $DATA_ROOT/02/deposit/sanntisv1_nucleotide_testsets_v3.tar.gz
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
here=$REPO/objectives/02-mibig31-validation
archives=$DATA_ROOT/01/inputs/mibig_5.0rc1-train
background=$DATA_ROOT/02/background
out=$DATA_ROOT/02/testsets
pool=$DATA_ROOT/01/pool
deposit=$DATA_ROOT/02/deposit/sanntisv1_nucleotide_testsets_v3.tar.gz

testsets() {
    [[ -s $archives/manifest.tsv ]] || die "no objective-01 archives in $archives (run objectives/01-production-models/01_inputs.sh)"
    [[ -s $background/background/manifest.json ]] || die "no background pool in $background (run 03_background.sh)"
    [[ -s $pool/v2/manifest.json ]] || die "no negative pool in $pool (run objectives/01-production-models/03_pool.sh)"
    "$PYTHON" "$here/workers/build_testsets.py" --archives "$archives" --pool "$pool" --pool-version v2 \
        --background "$background" --out "$out" \
        ${LIMIT:+--limit "$LIMIT"}
    local c=$out/contigs.tsv
    if [[ -z ${LIMIT:-} ]]; then
        expect_count "test BGCs (Test B)" "$("$PYTHON" -c 'import pandas as pd, sys; print(len(pd.read_parquet(sys.argv[1])))' "$out/test_bgcs.parquet")" 375
        expect_count "Test A BGCs" "$("$PYTHON" -c 'import pandas as pd, sys; print(int(pd.read_parquet(sys.argv[1])["test_a"].sum()))' "$out/test_bgcs.parquet")" 178
        expect_count "isolate contigs" "$(awk -F'\t' '$1 == "isolate"' "$c" | wc -l)" 70
        expect_count "shuffled contigs" "$(awk -F'\t' '$1 == "shuffled"' "$c" | wc -l)" 375
        expect_count "BGCs placed in the isolate set" "$(awk -F'\t' '$1 == "isolate" {s += $10} END {print s}' "$c")" 375
    fi
    mkdir -p "$(dirname "$deposit")"
    "$PYTHON" - "$out/testsets" "$deposit" <<'PY'
import sys
from pathlib import Path
from sv1.archives import sha256, write_tar
src, dest = Path(sys.argv[1]), Path(sys.argv[2])
write_tar(dest, {f"nucleotide_testsets/{p.relative_to(src)}": p.read_bytes() for p in src.rglob("*") if p.is_file()})
print(f"{sha256(dest)}  {dest.name}")
PY
}

step 02_testsets 2 16G 6:00:00 testsets
