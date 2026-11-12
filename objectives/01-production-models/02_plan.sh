#!/usr/bin/env bash
# The genome plan of negative pool v2: `sanntis train --until annotate` on 01_inputs.sh's archives (the production
# positives, deduplicated, with InterProScan 6 tokens), then `sanntis plan-negatives`: per-genome CDS quotas that
# follow the lineages of the positives (NCBI taxdump 2026-07-01 of objective 00), over the candidates of
# objective 00 (pool v1 genomes, BacDive type strains, RefSeq fungi) and the GTDB r226 species representatives,
# without the genomes of accessions/pool_exclude.tsv. The plan must equal objective 00's
# accessions/pool_v2_genomes.tsv, the table the pool genomes are downloaded by.
#   → $DATA_ROOT/01/plan/{run/, genomes.tsv, genomes.stats.json}
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need sanntis nextflow
d=$DATA_ROOT/01/inputs/mibig_5.0rc1-train
acc=$REPO/objectives/00-data-collection/accessions
out=$DATA_ROOT/01/plan

# Plan settings: negative CDS budgets per track (prokaryotes; fungi = 41.12 x the fungal positives' CDS), at most
# 1,500 CDS per genome, 250 CDS for one genome of each phylum without positives, quotas >= 100 CDS, seed 0; every
# pool v1 prokaryote not planned joins with <= 500 CDS.
PLAN_ARGS=(--budget proka=469554 --budget euk=186849 --cap 1500 --floor-cds 250 --min-quota 100 --seed 0
           --breadth-cds 500)

plan() {
    [[ -s $d/mibig_5.0rc1-train_gbk.tar.gz ]] || die "no training archives in $d (run 01_inputs.sh)"
    local tax
    tax=$DATA_ROOT/$(awk -F'\t' '$1 == "taxdump" {print $4}' "$acc/taxonomy.tsv")
    [[ -s $tax ]] || die "no taxdump at $tax (run objectives/00-data-collection/02_taxonomy.sh)"
    mapfile -t nf < <(nf_args)
    mapfile -t pos < <(train_archives "$d" mibig_5.0rc1-train)
    SANNTIS_THREADS=${SANNTIS_THREADS:-4} sanntis -v train --run "$out/run" --mibig none \
        --external-gbk "${pos[0]}" --external-json "${pos[1]}" \
        --container-engine "$CONTAINER_ENGINE" "${nf[@]}" --until annotate
    sanntis -v plan-negatives --out "$out/genomes.tsv" --positives "$out/run" --taxdump "$tax" \
        --candidates "$acc/pool_candidates_pool_v1.tsv" --candidates "$acc/pool_candidates_bacdive.tsv" \
        --candidates "$acc/pool_candidates_refseq_fungi.tsv" --gtdb-release 226 \
        --exclude "$acc/pool_exclude.tsv" --breadth "$acc/pool_breadth.tsv" "${PLAN_ARGS[@]}"
    # the plan of record (objective 00 downloads its genomes); a trial plans for its subset of positives
    [[ -z ${LIMIT:-} ]] || { log "LIMIT: plan not compared with accessions/pool_v2_genomes.tsv"; return 0; }
    diff <(cut -f1,2,5,10,18 "$acc/pool_v2_genomes.tsv") <(cut -f1,2,5,10,18 "$out/genomes.tsv") \
        || die "$out/genomes.tsv differs from accessions/pool_v2_genomes.tsv (genome, kingdom, accession, quota_cds, plan_role)"
    log "plan = accessions/pool_v2_genomes.tsv ✓"
}

step 01_plan 4 48G 2-00:00:00 plan
