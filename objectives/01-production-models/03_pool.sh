#!/usr/bin/env bash
# Negative pool v2: `sanntis build-negatives` on the genomes of the plan (objective 00 downloads), with the test
# backgrounds and held-out evaluation species of accessions/pool_exclude.tsv excluded: pyrodigal gene calls
# (prokaryotes) or RefSeq CDS (fungi), antiSMASH 8.0.1 `--minimal`, BGC-free runs, the isolate test backbones (55
# prokaryotic + 15 fungal runs, held out), per-genome portions of the BGC-free CDS (quota_cds) and InterProScan 6 on
# the portion and backbone proteins. Fan-out units (labels cpu_small, antismash) and InterProScan 6 go through
# Nextflow with NF_CONFIG; with RUNNER=slurm this step is the orchestrator job. Then the Hub deposit copy.
#   ANTISMASH_DATABASES=<dir>   antiSMASH 8.0.1 databases (the pinned image lacks some, e.g. transATor)
#   → $DATA_ROOT/01/pool/{genomes.tsv, prepared/, antismash/, ips6/, v2/}, $DATA_ROOT/01/hub/sanntis-training-data/negatives/v2/
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need sanntis nextflow
acc=$REPO/objectives/00-data-collection/accessions
asm=$DATA_ROOT/raw/ncbi/assemblies
out=$DATA_ROOT/01/pool

pool() {
    local table=$DATA_ROOT/01/pool_genomes.tsv
    # the plan with the downloaded genomes as local paths (prokaryotes: FASTA, fungi: annotated GenBank)
    "$PYTHON" - "$acc/pool_v2_genomes.tsv" "$asm" "$table" "${LIMIT:-0}" <<'PY'
import sys
from pathlib import Path
import pandas as pd
g = pd.read_csv(sys.argv[1], sep="\t", dtype=str, keep_default_na=False)
if int(sys.argv[4]):  # a trial: the first n genomes per kingdom (objective 00 downloads the same)
    g = g.groupby("kingdom", sort=False).head(int(sys.argv[4]))
ext = g["kingdom"].map({"proka": "fna.gz", "euk": "gbff.gz"})
g["path"] = [str(Path(sys.argv[2]) / a / f"{a}.{e}") for a, e in zip(g["accession"], ext, strict=True)]
missing = [p for p in g["path"] if not Path(p).is_file()]
if missing:
    sys.exit(f"{len(missing)} pool genomes not downloaded (objectives/00-data-collection/03_pool_genomes.sh), "
             f"e.g. {missing[0]}")
g.to_csv(sys.argv[3], sep="\t", index=False)
PY
    local args=(--container-engine "$CONTAINER_ENGINE") f
    for f in $NF_CONFIG; do args+=(--nf-config "$f"); done
    [[ -z ${ANTISMASH_DATABASES:-} ]] || args+=(--antismash-databases "$ANTISMASH_DATABASES")
    sanntis -v build-negatives --out "$out" --pool-version v2 --genomes "$table" --exclude "$acc/pool_exclude.tsv" \
        "${args[@]}"
    "$PYTHON" -c 'import json, sys
m = json.load(open(sys.argv[1]))
print(json.dumps({"counts": m["counts"], "excluded_genomes": m.get("excluded_genomes"),
                  "matches_api": m["ips"].get("matches_api")}, indent=1))' "$out/v2/manifest.json"
    [[ -n ${LIMIT:-} ]] || expect_count "isolate test backbones" "$("$PYTHON" -c 'import pandas as pd, sys
print(len(pd.read_parquet(sys.argv[1])))' "$out/v2/test_backbones.parquet")" 70
    # the Hub deposit (dataset microbiome-informatics/sanntis-training-data, negatives/v2/, tag negatives-v2)
    local hub=$DATA_ROOT/01/hub/sanntis-training-data/negatives/v2
    rm -rf "$hub" && mkdir -p "$(dirname "$hub")" && cp -r "$out/v2" "$hub"
    "$PYTHON" -c 'from sanntis.negpool.pool import load_pool; import sys; print(load_pool(sys.argv[1]).path, "sha256 ✓")' "$hub"
}

step 01_pool 4 64G 7-00:00:00 pool
