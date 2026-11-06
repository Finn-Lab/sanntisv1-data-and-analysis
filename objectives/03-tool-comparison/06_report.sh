#!/usr/bin/env bash
# Report: markdown + figures (metric tables with CIs, threshold sweep with the tools as points, strata radars,
# classification table).
#   → $DATA_ROOT/03/report/{report.md,figures/}
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
o=$REPO/objectives/03-tool-comparison

report() {
    "$PYTHON" "$o/workers/report.py" --score "$DATA_ROOT/03/score" --classification "$DATA_ROOT/03/classification" \
        --out "$DATA_ROOT/03/report"
}

step 03_report 1 4G 1:00:00 report
