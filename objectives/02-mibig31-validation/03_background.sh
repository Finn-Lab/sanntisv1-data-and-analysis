#!/usr/bin/env bash
# Test-background pool of the shuffled test set: objective 00's 124 test-background genomes (none in negative pool v2)
# through `sanntis build-negatives`' stages up to the antiSMASH stripping (gene calls, antiSMASH 8.0.1, native mask,
# BGC-free runs); no InterProScan. Fan-out units (labels cpu_small, antismash) go through
# Nextflow with NF_CONFIG and CONTAINER_ENGINE; with RUNNER=slurm this step is the orchestrator job.
#   ANTISMASH_DATABASES=<dir>   antiSMASH 8.0.1 databases; the pinned image lacks some (e.g. transATor). On Codon:
#                               /hps/nobackup/rdf/metagenomics/service-team/ref-dbs/antismash/8.0.1
#   BACKGROUND_ARGS="--runner inline --workers 8"   all units in this step's job instead (no Nextflow)
#   → $DATA_ROOT/02/background/background/{genomes,cds,clean_runs,native_mask,test_backbones}.parquet, manifest.json
#     (+ the build stages: prepared/<genome>/input.gbk, the backbone sequences and CDS of 04_testsets.sh)
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
[[ ${BACKGROUND_ARGS:-} == *"--runner inline"* ]] || need nextflow
here=$REPO/objectives/02-mibig31-validation
list=$REPO/objectives/00-data-collection/accessions/testset_background_genomes.tsv
out=$DATA_ROOT/02/background

background() {
    local args=(--container-engine "$CONTAINER_ENGINE") f
    for f in $NF_CONFIG; do args+=(--nf-config "$f"); done
    [[ -z ${ANTISMASH_DATABASES:-} ]] || args+=(--antismash-databases "$ANTISMASH_DATABASES")
    # shellcheck disable=SC2086  # BACKGROUND_ARGS: options, split on purpose
    local genomes=$DATA_ROOT/02/background_genomes.tsv
    { head -1 "$list"; rows_per_kingdom "$list" 3; } > "$genomes"  # LIMIT=n: the first n per kingdom
    "$PYTHON" "$here/workers/background_pool.py" "$genomes" "$DATA_ROOT/raw/ncbi/assemblies" "$out" "${args[@]}" \
        ${BACKGROUND_ARGS:-}
    # counts per kingdom (genomes, CDS, BGC-free runs, native CDS, backbones) = expected/
    "$PYTHON" -c 'import json, sys
c = json.load(open(sys.argv[1]))["counts"]
cols = ("genomes", "cds", "clean_runs", "native_cds", "test_backbones")
print("kingdom\t" + "\t".join(cols))
for k in ("proka", "euk"):
    print("\t".join([k] + [str(c[k][x]) for x in cols]))' "$out/background/manifest.json" > "$out/background/counts.tsv"
    column -t "$out/background/counts.tsv"
    if [[ -z ${LIMIT:-} && -s $here/expected/background_counts.tsv ]]; then
        diff "$here/expected/background_counts.tsv" "$out/background/counts.tsv" \
            || die "$out/background/counts.tsv differs from expected/background_counts.tsv"
    fi
}

step 02_background 2 16G 3-00:00:00 background
