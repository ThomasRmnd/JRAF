#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

#==============================
#  Utility logging function
#==============================

RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
BLUE='\033[0;34m'
NC='\033[0m'

log() {
    local level="$1"; shift
    local msg="$*"
    local timestamp
    timestamp="$(date '+%Y-%m-%d %H:%M:%S')"

    local level_num=0 color="$NC"
    case "$level" in
        ERROR) level_num=1; color="$RED" ;;
        WARN)  level_num=2; color="$YELLOW" ;;
        INFO)  level_num=3; color="$GREEN" ;;
        DEBUG) level_num=4; color="$CYAN" ;;
        ALL)   level_num=5; color="$BLUE" ;;
        *)     level_num=3 ;;
    esac

    local prefix="${color}[$timestamp][$level]${NC}"

    case "$level" in
        DEBUG|INFO) echo -e "${prefix} $msg" >&1 ;;
        WARN|ERROR) echo -e "${prefix} $msg" >&2 ;;
        ALL)
            echo -e "${prefix} $msg" >&1
            echo -e "${prefix} $msg" >&2
            ;;
        *) echo -e "${prefix} $msg" >&1 ;;
    esac
}

#==============================
# Utility functions
#==============================

HOSTNAME=$(hostname -f 2>/dev/null || hostname)
if [[ "${HOSTNAME}" =~ ^cc.*\.in2p3\.fr$ ]]; then
    # Detect CC-IN2P3 cluster
    CLUSTER="CC-IN2P3"
elif [[ "${HOSTNAME}" =~ ^lxlogin[0-9]+\.ihep\.ac\.cn$ ]]; then
    # Detect IHEP cluster
    CLUSTER="IHEP"
else
    echo "ERROR: Unknown cluster. Hostname: ${HOSTNAME}" >&2
    echo "Expected CC-IN2P3 (cca###) or IHEP (lxlogin###.ihep.ac.cn)" >&2
    exit 1
fi

log INFO "Cluster detected: ${CLUSTER}"

#==============================
# Configuration defaults
#==============================

usage() {
    cat <<EOF
Usage: $(basename "$0") --directory <path> --run <int>

Required:
    --directory <path>      Base directory
    --run       <int>       Run number

Optional:
    --help                  Show this help message and exit
EOF
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --directory) DIRECTORY="$2"; shift 2 ;;
            --run)       RUN_NUMBER="$2"; shift 2 ;;
            --help|-h)    usage; exit 0 ;;
            *) log ERROR "Unknown argument: $1"; usage; exit 1 ;;
        esac
    done

    if [[ -z "${DIRECTORY:-}" ]]; then
        log ERROR "--directory is required"
        usage
        exit 1
    fi

    if [[ -z "${RUN_NUMBER:-}" ]]; then
        log ERROR "--run is required"
        usage
        exit 1
    fi
}

#==============================
# Perform hadd for each run
#==============================

do_hadd() {
    local subdirectory="$1"
    local run="$2"
    local suffix="$3"

    local bucket_val=$(( (10#$run / 1000) * 1000 ))
    local group_val=$(( (10#$run / 100) * 100 ))
    
    local run_bucket=$(printf "%08d" "$bucket_val")
    local run_group=$(printf "%08d" "$group_val")

    local base_dir="${subdirectory}/${run_bucket}"

    mapfile -t group_candidates < <(find "${base_dir}" -maxdepth 1 -mindepth 1 -type d | grep -E "${run_group}(_[^/]+)?$" | sort)

    if (( ${#group_candidates[@]} == 0 )); then
        log WARN "No matching run_group directory found under ${base_dir} for group ${run_group}"
        return
    fi

    local selected_group_dir="${group_candidates[-1]}"

    local run_dir="${selected_group_dir}/${run}"

    if [[ ! -d "${run_dir}" ]]; then
        log WARN "Run directory not found: ${run_dir}"
        return
    fi

    local files=($(printf "%s\n" "${run_dir}"/*.root | sort -V))
    local nfiles=${#files[@]}
    if (( nfiles == 0 )); then
        log WARN "No ROOT files to merge for run ${run}"
        return
    fi


    local output_dir="${subdirectory}/summary"
    mkdir -p "${output_dir}"
    local output_file="${output_dir}/RUN.${run}.${suffix}"

    log INFO "Merging ${nfiles} files in ${run_dir} ..."
    hadd -f "${output_file}" "${files[@]}" && \
        log INFO "Merged successfully -> ${output_file}" || \
        log ERROR "hadd failed for run ${run}"
}

#==============================
# Main
#==============================

main() {
    parse_args "$@"
    do_hadd "${DIRECTORY}/analysis/ibd" "${RUN_NUMBER}" "analysis.root"
    do_hadd "${DIRECTORY}/reconstruction/reprod" "${RUN_NUMBER}" "reconstruction.root"
}

main "$@"