#!/usr/bin/env bash
# Report (markdown + figures: training, OOF validation, Test A / B metrics with CIs, sweeps, strata radars, NPC) and
# the v3active_eval deposit objective 03 reads (model bundle, test-BGC table, NPC labels, tokens, scores).
#   → $DATA_ROOT/02/report/{report.md,summary.tsv,figures/}, $DATA_ROOT/02/deposit/sanntisv1_v3active_eval_v3.tar.gz
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
here=$REPO/objectives/02-mibig31-validation
eval=$DATA_ROOT/02/eval

report() {
    [[ -s $eval/summary.json ]] || die "no scores at $eval (run 05_score.sh)"
    "$PYTHON" "$here/workers/report_v3active.py" --run "$DATA_ROOT/02/run" --eval "$eval" --out "$DATA_ROOT/02/report"
    # training side of the frozen tables (expected/; tune.tsv / validation_oof.tsv are derived from tune.json)
    compare_expected "$DATA_ROOT/02/run/bundle/report" "$here/expected" npc_validation.tsv sweep.tsv \
        validation_class_curves.tsv validation_strata.tsv
    expect_count "figures" "$(find "$DATA_ROOT/02/report/figures" -name '*.png' | wc -l)" 6
}

deposit_eval() {  # deterministic tarball: v3active_eval/{bundle/, test_bgcs.parquet, labels/, tokens.parquet, score/, summary.json}
    local out=$DATA_ROOT/02/deposit/sanntisv1_v3active_eval_v3.tar.gz
    mkdir -p "$DATA_ROOT/02/deposit"
    "$PYTHON" - "$eval" "$out" <<'PY'
import sys
from pathlib import Path
from sv1.archives import sha256, write_tar
src, dest = Path(sys.argv[1]), Path(sys.argv[2])
bundle = (src / "bundle").resolve()
members = {f"v3active_eval/bundle/{p.relative_to(bundle)}": p.read_bytes() for p in bundle.rglob("*") if p.is_file()}
for d in ("labels", "score"):
    members |= {f"v3active_eval/{p.relative_to(src)}": p.read_bytes() for p in (src / d).rglob("*") if p.is_file()}
for f in ("test_bgcs.parquet", "tokens.parquet", "summary.json"):
    members[f"v3active_eval/{f}"] = (src / f).read_bytes()
write_tar(dest, members)
print(f"{sha256(dest)}  {dest.name}  ({len(members)} members)")
PY
}

step 02_report 1 4G 1:00:00 report
step 02_deposit 1 4G 1:00:00 deposit_eval
