#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

#==============================
# Utility functions
#==============================

HOSTNAME=$(hostname -f 2>/dev/null || hostname)
if [[ "${HOSTNAME}" =~ ^cca[0-9]+\.in2p3\.fr$ ]]; then
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
# Configuration defaults
#==============================

RUN_LIST_REPROD25C="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/other/GoodList/ReProd25C/physics_good.txt"
RUN_LIST_REPROD25D="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/other/GoodList/ReProd25D/physics_good.txt"
RUN_LIST_REPROD26B="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/other/GoodList/ReProd26B/physics_good.txt"

OUTPUT_DIR="/sps/juno/jdeandre/rtraw_ThomasRaymond"

LOWER_BOUND=""
UPPER_BOUND=""

usage() {
    cat <<EOF
Usage: $(basename "$0") --campaign <str> [options]

Required:
  --campaign    <str>       Campaign selection {ReProd25C|ReProd25D|ReProd26B}

Optional:
  --lower       <num>       Starting run number (inclusive)
  --upper       <num>       Ending run number (inclusive)
  --output      <path>      Output path (default: ${OUTPUT_DIR})
  --help                    Show this help message and exit
EOF
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --campaign) CAMPAIGN="$2"; shift 2 ;;
            --lower)    LOWER_BOUND="$2"; shift 2 ;;
            --upper)    UPPER_BOUND="$2"; shift 2 ;;
            --output)   OUTPUT_DIR="$2"; shift 2 ;;
            --help|-h) usage; exit 0 ;;
            *) log ERROR "Unknown argument: $1"; usage; exit 1 ;;
        esac
    done

    if [[ -z "${CAMPAIGN:-}" ]]; then
        log ERROR "--campaign is required"
        usage
        exit 1
    fi

    case "${CAMPAIGN}" in
        ReProd25C)
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            ;;
        ReProd25D)
            RUN_LIST_PATH="${RUN_LIST_REPROD25D}"
            ;;
        ReProd26B)
            RUN_LIST_PATH="${RUN_LIST_REPROD26B}"
            ;;
        *)
            log ERROR "Invalid --campaign: ${CAMPAIGN}"
            usage
            exit 1
            ;;
    esac
}

#==============================
# Fetch Run List
#==============================

load_run_list() {
    log INFO "Fetching run list from EOS..."
    mapfile -t RUN_LIST < <(cat "${RUN_LIST_PATH}" | tr -d '\r' | sed '/^$/d')
}

#==============================
# Range Filtering
#==============================

filter_runs() {
    local filtered=()

    for run in "${RUN_LIST[@]}"; do
        (( run < 0 )) && continue
        if [[ -n "${LOWER_BOUND}" && "${run}" -lt "${LOWER_BOUND}" ]]; then
            continue
        fi
        if [[ -n "${UPPER_BOUND}" && "${run}" -gt "${UPPER_BOUND}" ]]; then
            continue
        fi
        filtered+=("${run}")
    done

    RUN_LIST=("${filtered[@]}")

    if (( ${#RUN_LIST[@]} == 0 )); then
        log WARN "No runs matched the provided range"
        exit 0
    fi

    log INFO "Runs selected: ${RUN_LIST[*]}"
}

#==============================
# Main
#==============================

main() {
    parse_args "$@"
    load_run_list
    filter_runs

    for run in "${RUN_LIST[@]}"; do
        log INFO "Submitting job for run ${run}"

        sbatch \
            --job-name="hadd_${run}_batch" \
            --output="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/log/hadd_${run}.log" \
            --error="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/err/hadd_${run}.err" \
            --partition="htc" \
            --ntasks=1 \
            --cpus-per-task=1 \
            --mem="1G" \
            --time="0-0:02:30" \
            --mail-user="thomas.raymond@iphc.cnrs.fr" \
            --mail-type="FAIL" \
            job_worker.sh \
            --directory "${OUTPUT_DIR}" --run "${run}"
    done

    log INFO "All jobs submitted successfully"
}

main "$@"