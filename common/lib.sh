# shellcheck shell=bash
# Shared helpers for the runbooks. Sourced by common/env.sh; nothing here runs on its own.

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

need() {  # need <command>... : fail early when a tool is missing
    local c
    for c in "$@"; do command -v "$c" >/dev/null || die "'$c' not found on PATH"; done
}

# ── downloads (every input is fetched from its repository and checked) ────────

fetch() {  # fetch <url> <dest> [sha256] : download once; verify when a checksum is given
    local url=$1 dest=$2 sha=${3:-}
    if [[ ! -s $dest ]]; then
        mkdir -p "$(dirname "$dest")"
        log "fetch ${url%%&api_key=*}"  # never log the NCBI key
        curl -fsSL --retry 5 --retry-delay 10 -o "$dest.part" "$url" || die "download failed: ${url%%&api_key=*}"
        mv "$dest.part" "$dest"
    fi
    [[ -z $sha ]] || echo "$sha  $dest" | sha256sum -c --quiet - \
        || { rm -f "$dest"; die "sha256 mismatch, file removed (re-run to download it again): $dest"; }
}

fetch_zenodo() {  # fetch_zenodo <record id> <file name> <dest dir> <sha256>  (ZENODO_MIRROR=<dir>: local copies)
    if [[ -n ${ZENODO_MIRROR:-} ]]; then
        local src
        src=$(find "$ZENODO_MIRROR" -name "$2" -print -quit)
        [[ -n $src ]] || die "$2 not in ZENODO_MIRROR ($ZENODO_MIRROR)"
        mkdir -p "$3"
        echo "$4  $3/$2" | sha256sum -c --quiet - >/dev/null 2>&1 || cp -f "$src" "$3/$2"
        echo "$4  $3/$2" | sha256sum -c --quiet - || { rm -f "$3/$2"; die "sha256 mismatch: $src"; }
        return 0
    fi
    [[ $1 != TBD ]] || die "$2: Zenodo record not published yet (set ZENODO_MIRROR to a local copy)"
    fetch "https://zenodo.org/records/$1/files/$2?download=1" "$3/$2" "$4"
}

deposit() {  # deposit <name> <dest dir> : fetch + unpack a Zenodo deposit listed in env/zenodo.tsv
    local rec file sha
    read -r rec file sha < <(awk -F'\t' -v d="$1" '$1 == d {print $2, $3, $4; f=1} END {exit !f}' "$REPO/env/zenodo.tsv") \
        || die "deposit '$1' not in env/zenodo.tsv"
    [[ -e $2/.unpacked-$sha ]] && return 0
    fetch_zenodo "$rec" "$file" "$DATA_ROOT/zenodo" "$sha"
    mkdir -p "$2" && tar -xzf "$DATA_ROOT/zenodo/$file" -C "$2" && touch "$2/.unpacked-$sha"
}

fetch_assembly() {  # fetch_assembly <GCA_/GCF_ accession> <dest dir> [prot|gbff] : <dest>/<acc>/<acc>.fna.gz (+ .gff.gz)
    # via NCBI Datasets; with "prot" also <acc>.faa.gz, the annotated proteins, when the assembly has an annotation;
    # with "gbff" also <acc>.gbff.gz, the annotated GenBank (RefSeq CDS of eukaryotic backbones)
    local acc=$1 out=$2/$1 types=GENOME_FASTA,GENOME_GFF
    [[ ${3:-} == prot ]] && types+=,PROT_FASTA
    [[ ${3:-} == gbff ]] && types+=,GENOME_GBFF
    [[ -e $out/done ]] && return 0
    rm -rf "$out" && mkdir -p "$out"
    fetch "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=$types" \
        "$out/$acc.zip"
    unzip -qo "$out/$acc.zip" -d "$out" && rm "$out/$acc.zip"
    (cd "$out" && md5sum -c --quiet md5sum.txt) || die "md5 mismatch in the NCBI Datasets package of $acc"
    local d=$out/ncbi_dataset/data
    ls "$d/$acc"/*_genomic.fna >/dev/null 2>&1 || die "no genome FASTA for $acc"
    gzip -c "$d/$acc"/*_genomic.fna > "$out/$acc.fna.gz"
    [[ ! -e $d/$acc/genomic.gff ]] || gzip -c "$d/$acc/genomic.gff" > "$out/$acc.gff.gz"
    [[ ! -e $d/$acc/protein.faa ]] || gzip -c "$d/$acc/protein.faa" > "$out/$acc.faa.gz"
    [[ ! -e $d/$acc/genomic.gbff ]] || gzip -c "$d/$acc/genomic.gbff" > "$out/$acc.gbff.gz"
    [[ ${3:-} != gbff || -s $out/$acc.gbff.gz ]] || die "no annotated GenBank (GBFF) for $acc"
    mv "$d/assembly_data_report.jsonl" "$out/"
    rm -rf "$out/ncbi_dataset" "$out/README.md" "$out/md5sum.txt"
    touch "$out/done"
}

each() {  # each <jobs> <function> : call <function> with the fields of every stdin line, <jobs> at a time
    local jobs=$1 fn=$2
    # shellcheck disable=SC2046  # every function (helpers and the runbook's) is visible to the workers
    export -f $(declare -F | awk '{print $3}')
    # shellcheck disable=SC2016  # $0 / $@ expand in the child bash
    xargs -P "$jobs" -L 1 bash -c 'set -euo pipefail; "$0" "$@"' "$fn"
}

sums() {  # sums <dir> <glob> : write <dir>/SHA256SUMS for the files matching <glob> (record of what was fetched)
    (cd "$1" && find . -type f -name "$2" -print0 | sort -z | xargs -0 -r sha256sum > SHA256SUMS)
}

rows() {  # rows <table> : data rows of a TSV (LIMIT=n keeps the first n: a trial run, use a separate DATA_ROOT)
    awk -v n="${LIMIT:-0}" 'NR > 1 && (n == 0 || NR <= n + 1)' "$1"  # awk, not head: no SIGPIPE under pipefail
}

rows_per_kingdom() {  # rows_per_kingdom <table> <kingdom column> : like rows, LIMIT=n keeps the first n per kingdom
    awk -F'\t' -v n="${LIMIT:-0}" -v k="$2" 'NR > 1 && (n == 0 || ++c[$k] <= n)' "$1"
}

train_archives() {  # train_archives <dir> <name> : the <name>_{gbk,json}.tar.gz pair of <dir> (LIMIT=n: a subset of
    # 10 x n BGCs, written to <dir>/limit/), one path per line
    local g=$1/$2_gbk.tar.gz j=$1/$2_json.tar.gz
    if [[ -n ${LIMIT:-} ]]; then
        # + the first 2 x n entries >= BGC0002750: they hold a trial's first n test BGCs (NPC labels in collect)
        "$PYTHON" -m sv1.archives limit "$g" "$j" "$1/limit" "$((LIMIT * 10))" "$((LIMIT * 2))"
    else
        printf '%s\n%s\n' "$g" "$j"
    fi
}

compare_expected() {  # compare_expected <new dir> <expected dir> <file>... : the frozen tables of a full run
    # (expected/) vs this run's, with sv1.expected (6 significant digits); not in a LIMIT trial, nor before freezing
    [[ -z ${LIMIT:-} ]] || { log "LIMIT: $2 not compared"; return 0; }
    local f missing=()
    for f in "${@:3}"; do [[ -e $2/$f ]] || missing+=("$f"); done
    if ((${#missing[@]})); then log "not frozen in $2 (no comparison): ${missing[*]}"; return 0; fi
    "$PYTHON" -m sv1.expected compare "$1" "$2" --files "${@:3}"
}

expect_count() {  # expect_count <label> <found> <expected>
    [[ $2 -eq $3 ]] || die "$1: $2 found, $3 expected"
    log "$1: $2 / $3 ✓"
}

fetch_nuccore() {  # fetch_nuccore <accession> <from> <to> <dest gbk> : one locus as GenBank via NCBI E-utilities
    local q="db=nuccore&id=$1&rettype=gbwithparts&retmode=text"
    [[ $3 -gt 0 ]] && q+="&seq_start=$2&seq_stop=$3"  # 0/0 in MIBiG = the whole record
    [[ -s $4 ]] && return 0
    fetch "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?$q${NCBI_API_KEY:+&api_key=$NCBI_API_KEY}" "$4"
    grep -q '^//' "$4" || { rm -f "$4"; die "truncated GenBank from nuccore: $1"; }
    sleep "${NCBI_SLEEP:-0.4}"  # E-utilities rate limit (3/s without a key)
}

# ── containers ────────────────────────────────────────────────────────────────

image() {  # image <tool> : pinned image from env/containers*.tsv, or IMAGE_<tool> (e.g. IMAGE_gecco_0_9_5=<.sif>)
    # env/containers.tsv, then any env/containers.<name>.tsv a derived repository adds for its own tools
    local v="IMAGE_${1//[^A-Za-z0-9]/_}" img
    if [[ -n ${!v:-} ]]; then echo "${!v}"; return 0; fi
    img=$(awk -F'\t' -v t="$1" '$1 == t {print $3; f=1; exit} END {exit !f}' "$REPO"/env/containers*.tsv) \
        || die "no image for '$1' in env/containers*.tsv"
    [[ $img != TODO* ]] || die "$1: image not published yet (${img#TODO: }); build it and set $v"
    echo "$img"
}

ctr() {  # ctr <tool> <command>... : run a pinned tool container with DATA_ROOT, the cwd and CTR_BINDS mounted
    local img b binds=("$DATA_ROOT" "$PWD") args=()
    img=$(image "$1"); shift
    if [[ $img == host ]]; then "$@"; return; fi  # IMAGE_<tool>=host: the tool is installed on this machine
    for b in ${CTR_BINDS:-}; do binds+=("$b"); done  # extra mounts: <dir> or <host dir>:<container dir>
    case $CONTAINER_ENGINE in
        docker)
            for b in "${binds[@]}"; do [[ $b == *:* ]] && args+=(-v "$b") || args+=(-v "$b:$b"); done
            [[ -z ${TORCH_HOME:-} ]] || args+=(-e "TORCH_HOME=$TORCH_HOME")
            docker run --rm -u "$(id -u):$(id -g)" "${args[@]}" -w "$PWD" "$img" "$@" ;;
        singularity|apptainer)  # a local image (.sif / sandbox) or a registry reference; the environment passes
            for b in "${binds[@]}"; do args+=(-B "$b"); done
            [[ -e $img ]] || img=docker://$img
            "$CONTAINER_ENGINE" exec "${args[@]}" "$img" "$@" ;;
        *) die "CONTAINER_ENGINE must be docker, singularity or apptainer (got '$CONTAINER_ENGINE')" ;;
    esac
}

nf_args() {  # nf_args : --nf-config <file> for each file in NF_CONFIG
    local f
    for f in $NF_CONFIG; do printf -- '--nf-config\n%s\n' "$f"; done
}

# ── steps ─────────────────────────────────────────────────────────────────────

step() {  # step <name> <cpus> <mem> <time> <command>... : run once (done marker), locally or as a SLURM job
    local name=$1 cpus=$2 mem=$3 time=$4; shift 4
    local mark=$WORK/.done/$name log=$WORK/logs/$name.log
    if [[ -n ${STEP_ONLY:-} ]]; then  # inside a SLURM job: run only the requested step, here
        [[ $name == "$STEP_ONLY" ]] || return 0
    elif [[ -e $mark && -z ${FORCE:-} ]]; then
        log "skip $name (done)"; return 0
    fi
    mkdir -p "$WORK/.done" "$WORK/logs"
    log "run $name"
    if [[ $RUNNER == local || -n ${STEP_ONLY:-} ]]; then
        # not `if ! cmd`: bash ignores errexit inside a condition, so a failing command would not stop the step
        local rc
        set +e
        (set -e; "$@") 2>&1 | tee "$log"
        rc=${PIPESTATUS[0]}
        set -e
        [[ $rc -eq 0 ]] || die "$name failed (log: $log)"
    elif [[ $RUNNER == slurm ]]; then
        # the job re-runs this runbook with the same arguments and STEP_ONLY=<name>, so it rebuilds the same context
        rm -f "$log"
        # shellcheck disable=SC2086
        sbatch --wait --parsable -J "$name" -c "$cpus" --mem "$mem" -t "$time" ${SLURM_ARGS:-} \
            -o "$WORK/logs/$name.slurm.log" \
            --wrap "STEP_ONLY=$name bash $(printf '%q ' "$SCRIPT" "${SCRIPT_ARGS[@]}")" >/dev/null \
            || die "$name failed (logs: $log, $WORK/logs/$name.slurm.log)"
        # a job that cannot write its log (WORK not on a shared filesystem) can still exit 0
        [[ -e $log ]] || die "$name: no job log; is \$WORK ($WORK) on a shared filesystem?"
    else
        die "RUNNER must be local or slurm (got '$RUNNER')"
    fi
    touch "$mark"
}
