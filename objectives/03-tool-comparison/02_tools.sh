#!/usr/bin/env bash
# Run the comparison tools on the nucleotide test sets, from the pinned containers (env/containers.tsv).
#   → $DATA_ROOT/03/tools/<run>/<set>_<track>/<contig>.<ext>   (the layout 04_score.sh reads)
#   → $DATA_ROOT/03/ips576/<set>/<track>/<contig>.ips576.tsv.gz   (run ips576: SanntiS 0.9.5's input)
# Usage: 02_tools.sh [run ...]   runs: gecco-0.9.5 gecco-0.11.0 ips576 sanntis-0.9.5 bgc-prophet antismash-relaxed
#                                      antismash-strict antismash-loose plantismash bgcat   (default: all)
# Every tool gets the same input: the contig GenBank with its CDS calls (gene finders that ignore CDS features
# call genes themselves). TOOL_CPUS per tool call, JOBS calls at a time; LIMIT=n: first n
# contigs per set and track (a trial, use a separate DATA_ROOT).
# shellcheck source=SCRIPTDIR/../../common/env.sh
source "$(dirname "$0")/../../common/env.sh"
need tar gzip
w=$REPO/objectives/03-tool-comparison/workers
in=$DATA_ROOT/03/inputs
tools=$DATA_ROOT/03/tools
SETS=(isolate shuffled)
export TOOL_CPUS=${TOOL_CPUS:-4} w in
export IPS_SH=${IPS_SH:-/opt/interproscan/interproscan.sh}  # the InterProScan 5.76 launcher inside the image (or host)
export IPS_TIMEOUT=${IPS_TIMEOUT:-2h}  # per contig attempt (a few hundred proteins: minutes)
export APPL=Pfam,NCBIfam,SFLD,SUPERFAMILY,Gene3D,Hamap,Coils,CDD,PRINTS,PIRSF,ProSiteProfiles,ProSitePatterns,MobiDBLite,SMART

scratch() { mktemp -d -p "$WORK/tmp"; }

# IMAGE_<tool>=host: <TOOL>_ENV = the prefix of the tool's environment (its bin/ goes first on PATH, for the tool
# and the executables it calls); <TOOL>_HOME = a source checkout where the image has one under /opt
host_env() { [[ -z $1 ]] || export PATH=$1/bin:$PATH; }

gecco_one() {  # <version> <gbk> <out dir>
    local p d
    p=$(basename "$2" .gbk); d=$(scratch)
    [[ -e $3/$p.done ]] && return 0
    (cd "$d" && ctr "gecco-$1" gecco run -j "$TOOL_CPUS" -o res -g "$2" > gecco.log 2>&1) \
        || die "GECCO $1 failed on $p ($d/gecco.log)"
    [[ ! -s $d/res/$p.clusters.tsv ]] || cp "$d/res/$p.clusters.tsv" "$3/"
    touch "$3/$p.done"; rm -rf "$d"
}

ips576_one() {  # <cds.parquet> <out tsv.gz> : InterProScan 5.76-107.0 TSV of the contig's proteins (ids = /protein_id)
    local d
    [[ -s $2 ]] && return 0
    d=$(mktemp -d -p "${TMPDIR:-$WORK/tmp}")  # node-local under SLURM: IPS temp files are many and small
    "$PYTHON" "$w/proteins.py" "$1" "$d/in.faa" cds_id
    for a in 1 2; do  # IPS rarely hangs in its output step; a re-run succeeds
        (cd "$d" && ctr interproscan-5.76 timeout "$IPS_TIMEOUT" "$IPS_SH" -i in.faa -f TSV -o out.tsv \
            -cpu "$TOOL_CPUS" -T "$d/t" --disable-precalc --iprlookup --goterms -appl "$APPL" > ips.log 2>&1) && break
        [[ $a -eq 2 ]] && { cp "$d/ips.log" "$2.ips.log"; rm -rf "$d"
                            die "InterProScan 5.76 failed twice on $(basename "$1" .cds.parquet) ($2.ips.log)"; }
        rm -rf "$d/t" "$d/out.tsv"
    done
    gzip -c "$d/out.tsv" > "$2.part"  # plain commands: errexit stops on a failed write (an && list would not)
    mv "$2.part" "$2"; rm -rf "$d"
}

sanntis095_one() {  # <gbk> <out dir> <ips 5.76 tsv.gz>
    local p d
    p=$(basename "$1" .gbk); d=$(scratch)
    [[ -e $2/$p.sanntis.gff ]] && return 0
    zcat "$3" > "$d/ips.tsv"
    (cd "$d" && ctr sanntis-0.9.5 sanntis --ip-file ips.tsv --outdir out --outfile "$p.sanntis.gff" "$1" > s.log 2>&1) \
        || die "SanntiS 0.9.5 failed on $p ($d/s.log)"
    cp "$d/$p.sanntis.gff" "$2/" 2>/dev/null || cp "$d/out/$p.sanntis.gff" "$2/"
    rm -rf "$d"
}

antismash_one() {  # <strictness> <gbk> <out dir> <taxon>
    local p d
    p=$(basename "$2" .gbk); d=$(scratch)
    [[ -s $3/$p.json.gz ]] && return 0
    "$PYTHON" "$w/gbk_ids.py" shorten "$2" "$d/in.gbk" "$d/ids.tsv"
    (cd "$d" && ctr antismash antismash --taxon "$4" -c "$TOOL_CPUS" ${ANTISMASH_DB:+--databases "$ANTISMASH_DB"} \
        --hmmdetection-strictness "$1" --genefinding-tool none --output-dir res --output-basename "$p" in.gbk \
        > antismash.log 2>&1) || die "antiSMASH $1 failed on $p ($d/antismash.log)"
    "$PYTHON" "$w/gbk_ids.py" restore-json "$d/res/$p.json" "$d/ids.tsv"
    gzip -c "$d/res/$p.json" > "$3/$p.json.gz"; rm -rf "$d"
}

plantismash_one() {  # <gbk> <out dir>
    local p d res
    host_env "${PLANTISMASH_ENV:-}"
    p=$(basename "$1" .gbk); d=$(scratch)
    [[ -s $2/$p.gbk.gz ]] && return 0
    "$PYTHON" "$w/gbk_ids.py" shorten "$1" "$d/in.gbk" "$d/ids.tsv"
    (cd "$d" && ctr plantismash python "${PLANTISMASH_HOME:-/opt/plantismash}/run_antismash.py" --taxon plants \
        -c "$TOOL_CPUS" --disable-html --disable-svg --disable-xls --genefinding none --outputfolder res in.gbk \
        > p.log 2>&1) \
        || die "plantiSMASH failed on $p ($d/p.log)"
    res=$(find "$d/res" -name '*.final.gbk' | head -1)
    [[ -n $res ]] || die "plantiSMASH: no final GenBank for $p"
    cp "$d/ids.tsv" "$2/$p.ids.tsv"; gzip -c "$res" > "$2/$p.gbk.gz"; rm -rf "$d"
}

contigs() {  # <set> <track> <out dir> : "<gbk> <out dir>" per contig of that set and track
    find "$in/nucleotide_testsets/$1/$2" -name '*.gbk' | sort | awk -v o="$3" -v n="${LIMIT:-0}" 'n == 0 || NR <= n {print $1, o}'
}

gecco() {  # <version>
    local s t o
    for s in "${SETS[@]}"; do for t in proka euk; do
        o=$tools/gecco_$1/${s}_$t; mkdir -p "$o"
        contigs "$s" "$t" "$o" | awk -v v="$1" '{print v, $0}' | each "${JOBS:-8}" gecco_one
    done; done
}

ips576() {  # InterProScan 5.76-107.0 on every contig's proteins (SanntiS 0.9.5's --ip-file input)
    local s t o
    for s in "${SETS[@]}"; do for t in proka euk; do
        o=$DATA_ROOT/03/ips576/$s/$t; mkdir -p "$o"
        contigs "$s" "$t" "$o" | while read -r g d; do
            echo "${g%.gbk}.cds.parquet $d/$(basename "$g" .gbk).ips576.tsv.gz"
        done | each "${JOBS:-8}" ips576_one
    done; done
}

sanntis095() {
    local s t o f
    for s in "${SETS[@]}"; do for t in proka euk; do
        o=$tools/sanntis_0.9.5.0/${s}_$t; mkdir -p "$o"
        contigs "$s" "$t" "$o" | while read -r g d; do
            f=$DATA_ROOT/03/ips576/$s/$t/$(basename "$g" .gbk).ips576.tsv.gz
            [[ -s $f ]] || die "no InterProScan 5.76 table for $(basename "$g" .gbk) (run 02_tools.sh ips576 first)"
            echo "$g $d $f"
        done | each "${JOBS:-8}" sanntis095_one
    done; done
}

antismash() {  # <strictness> : bacteria on prokaryotic contigs, fungi (fungiSMASH) on eukaryotic ones
    local s t o
    for s in "${SETS[@]}"; do for t in proka euk; do
        o=$tools/antismash_7.1.0_$1/${s}_$t; mkdir -p "$o"
        contigs "$s" "$t" "$o" | awk -v l="$1" -v x="$([[ $t == proka ]] && echo bacteria || echo fungi)" \
            '{print l, $1, $2, x}' | each "${JOBS:-8}" antismash_one
    done; done
}

plantismash() {  # eukaryotic contigs (plant BGCs sit on fungal / protist backbones)
    local s o
    for s in "${SETS[@]}"; do
        o=$tools/plantismash_1.0/${s}_euk; mkdir -p "$o"
        contigs "$s" euk "$o" | each "${JOBS:-8}" plantismash_one
    done
}

# Third-party weights are fetched from their publishers, never mirrored: BGC-Prophet's models (MIT; published only
# as a GitHub issue attachment) and the ESM2 8M weights fair-esm loads (MIT), pinned by sha256.
PROPHET_MODEL_URL=https://github.com/HUST-NingKang-Lab/BGC-Prophet/files/12733164/model.tar.gz
PROPHET_MODEL_SHA=bb4bff5e14c950f7c02884657d90536db9437787bd14d339f6dbc5fef3cffe08
ESM2_URL=https://dl.fbaipublicfiles.com/fair-esm
ESM2_SHA=46f002a9870c9bdecd0ea887acb1f9a38a6b561e8f8bf8a6990b679b9d31b928
ESM2_REG_SHA=8f7a4557d57713b97ba0e484303007efb7230d25299c0ac47a0a1b12a87bbb9d

prophet() {  # BGC-Prophet 0.1.2 detection on each contig's proteins in gene order
    local wt=$DATA_ROOT/cache/bgc_prophet ck=$DATA_ROOT/cache/torch/hub/checkpoints s t o faa d c
    fetch "$PROPHET_MODEL_URL" "$wt/model.tar.gz" "$PROPHET_MODEL_SHA"
    [[ -e $wt/model/annotator.pt ]] || tar -xzf "$wt/model.tar.gz" -C "$wt"
    fetch "$ESM2_URL/models/esm2_t6_8M_UR50D.pt" "$ck/esm2_t6_8M_UR50D.pt" "$ESM2_SHA"
    fetch "$ESM2_URL/regression/esm2_t6_8M_UR50D-contact-regression.pt" \
        "$ck/esm2_t6_8M_UR50D-contact-regression.pt" "$ESM2_REG_SHA"
    export TORCH_HOME=$DATA_ROOT/cache/torch
    host_env "${PROPHET_ENV:-}"
    for s in "${SETS[@]}"; do for t in proka euk; do
        o=$tools/bgc-prophet_0.1.2/${s}_$t; faa=$DATA_ROOT/03/prophet_in/${s}_$t
        [[ -s $o/prophet.csv ]] && continue
        mkdir -p "$o" "$faa"
        for c in "$in/nucleotide_testsets/$s/$t"/*.cds.parquet; do
            "$PYTHON" "$w/proteins.py" "$c" "$faa/$(basename "$c" .cds.parquet).faa"
        done
        d=$(scratch)
        # the `pipeline` subcommand's steps without its final classification, which reopens the LMDB of the
        # prediction step and fails with py-lmdb 3 ("already open in this process"); detection output identical
        (cd "$d" && ctr bgc-prophet bgc_prophet extract esm2_t6_8M_UR50D "$faa" "$d/lmdb" --directory --include mean \
                --toks_per_batch 40960 \
            && ctr bgc-prophet bgc_prophet organize --genomesDir "$faa" --outputPath "$d" --name genomes \
                --threads "$TOOL_CPUS" \
            && ctr bgc-prophet bgc_prophet split --genomesPath "$d/genomes.csv" --outputPath "$d" --name genomes \
                --threads "$TOOL_CPUS" \
            && ctr bgc-prophet bgc_prophet predict --datasetPath "$d/genomes_split.csv" --outputPath "$d" \
                --modelPath "$wt/model/annotator.pt" --lmdbPath "$d/lmdb" --name prophet \
                --device "${PROPHET_DEVICE:-cpu}" --saveIntermediate \
            && ctr bgc-prophet bgc_prophet output --datasetPath "$d/genomes_split.csv" --outputPath "$d" \
                --loadIntermediate "$d/intermediate_prophet.npy" --name prophet --threshold 0.5 --max_gap 3 \
                --min_count 2) > "$d/prophet.log" 2>&1 || die "BGC-Prophet failed on ${s}_$t ($d/prophet.log)"
        cp "$d/prophet.csv" "$o/prophet.csv"; rm -rf "$d"
    done; done
}

ESMC_URL=https://huggingface.co/biohub/esmc-600m-2024-12/resolve/e4d83bc7e10fd55c92e598e545f4a76bf04a6e5c/data/weights/esmc_600m_2024_12_v0.pth
ESMC_SHA=8ef856e1a237ee3f995442df997a962e70057faadecf38fc0c8561bd3c2f4324

bgcat() {  # BGCat on the biosynthetic genes antiSMASH 7.1.0 (relaxed) finds in each test BGC (isolate set)
    local as=$tools/antismash_7.1.0_relaxed o=$tools/bgcat_4f77517 wt=$DATA_ROOT/cache/esmc
    [[ -d $as/isolate_proka ]] || die "no antiSMASH 7.1.0 relaxed run under $as (run 02_tools.sh antismash-relaxed first)"
    mkdir -p "$o"
    "$PYTHON" "$w/bgcat_input.py" --testsets "$in/nucleotide_testsets" --antismash "$as" --out "$o/bgcs.json"
    fetch "$ESMC_URL" "$wt/esmc_600m_2024_12_v0.pth" "$ESMC_SHA"  # ESM C 600M (EvolutionaryScale licence)
    local home=${BGCAT_HOME:-/opt/bgcat}  # the repository inside the image (or on this machine with IMAGE_bgcat=host)
    host_env "${BGCAT_ENV:-}"            # IMAGE_bgcat=host: weights go in $home/data/weights (no mount)
    CTR_BINDS="${CTR_BINDS:-} $wt:$home/data/weights" ctr bgcat sh -c \
        "cd $home && python embed.py $o/bgcs.json $o/features.json && python classify.py $o/features.json $o/predictions.json" \
        > "$o/bgcat.log" 2>&1 || die "BGCat failed ($o/bgcat.log)"
}

mkdir -p "$WORK/tmp"
runs=("$@")
[[ ${#runs[@]} -gt 0 ]] || runs=(gecco-0.9.5 gecco-0.11.0 ips576 sanntis-0.9.5 bgc-prophet antismash-relaxed
                                 antismash-strict antismash-loose plantismash bgcat)
for r in "${runs[@]}"; do
    case $r in
        gecco-*)        step "03_tool_$r" 8 16G 12:00:00 gecco "${r#gecco-}" ;;
        ips576)         step "03_tool_$r" 32 128G 48:00:00 ips576 ;;  # JOBS=8 x TOOL_CPUS=4, ~14 GB JVM each
        sanntis-0.9.5)  step "03_tool_$r" 8 32G 12:00:00 sanntis095 ;;
        bgc-prophet)    step "03_tool_$r" 8 32G 12:00:00 prophet ;;
        antismash-*)    step "03_tool_$r" 16 64G 24:00:00 antismash "${r#antismash-}" ;;
        plantismash)    step "03_tool_$r" 8 32G 12:00:00 plantismash ;;
        bgcat)          step "03_tool_$r" 8 64G 24:00:00 bgcat ;;
        *) die "unknown run '$r'" ;;
    esac
done
