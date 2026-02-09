#!/bin/bash
#--------------------------------------------------------------------------------------------------
#  JUNO Multi-run Hadd Helper
#  Purpose: hadd all ROOT outputs for multiple runs automatically
#--------------------------------------------------------------------------------------------------

set -euo pipefail
IFS=$'\n\t'

#==============================
# Utility functions
#==============================

HOSTNAME=$(hostname -f 2>/dev/null || hostname)
if [[ "${HOSTNAME}" =~ ^cca[0-9]+\.in2p3\.fr$ ]]; then
    # Detect CC-IN2P3 cluster
    CLUSTER="CC-IN2P3"
    source /pbs/home/t/traymond/share/bash/logging.sh
elif [[ "${HOSTNAME}" =~ ^lxlogin[0-9]+\.ihep\.ac\.cn$ ]]; then
    # Detect IHEP cluster
    CLUSTER="IHEP"
    source /junofs/users/traymond/bash/logging.sh
else
    echo "ERROR: Unknown cluster. Hostname: ${HOSTNAME}" >&2
    echo "Expected CC-IN2P3 (cca###) or IHEP (lxlogin###.ihep.ac.cn)" >&2
    exit 1
fi

log INFO "Cluster detected: ${CLUSTER}"

#==============================
# Configuration defaults
#==============================

XRD_URL_EOS="root://junoeos01.ihep.ac.cn/"
RUN_LIST_REPROD25C="/eos/juno/groups/DataQuality/P25A/Physics/goodrunlist_v3.6/Physics_good_run_list.txt"
RUN_LIST_REPROD25D="/eos/juno/groups/DataQuality/ReProd25D/Physics/goodrunlist_v0.0-v2/physics_good.txt"

OUTPUT_DIR_BASE="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod"
OUTPUT_SUFFIX_NORMAL="output.normal.cca.root"
OUTPUT_SUFFIX_REPROD25A="output.reprod25a.cca.root"
OUTPUT_SUFFIX_REPROD25B="output.reprod25b.cca.root"
OUTPUT_SUFFIX_REPROD25C="output.reprod25c.cca.root"
OUTPUT_SUFFIX_REPROD25D="output.reprod25d.cca.root"

LOWER_BOUND=""
UPPER_BOUND=""

usage() {
    cat <<EOF
Usage: $(basename "$0") --campaign <str> [options]

Required:
  --campaign <str>      Campaign selection {Normal|ReProd25A|ReProd25B|ReProd25C|ReProd25D}

Optional:
  --lower <num>         Starting run number (inclusive)
  --upper <num>         Ending run number (inclusive)
  --help                Show this help message and exit
EOF
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --campaign) CAMPAIGN="$2"; shift 2 ;;
            --lower)    LOWER_BOUND="$2"; shift 2 ;;
            --upper)    UPPER_BOUND="$2"; shift 2 ;;
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
        Normal)
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            OUTPUT_SUFFIX="${OUTPUT_SUFFIX_NORMAL}"
            ;;
        ReProd25A)
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            OUTPUT_SUFFIX="${OUTPUT_SUFFIX_REPROD25A}"
            ;;
        ReProd25B)
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            OUTPUT_SUFFIX="${OUTPUT_SUFFIX_REPROD25B}"
            ;;
        ReProd25C)
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            OUTPUT_SUFFIX="${OUTPUT_SUFFIX_REPROD25C}"
            ;;
        ReProd25D)
            RUN_LIST_PATH="${RUN_LIST_REPROD25D}"
            OUTPUT_SUFFIX="${OUTPUT_SUFFIX_REPROD25D}"
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
    mapfile -t RUN_LIST < <(xrdfs "${XRD_URL_EOS}" cat "${RUN_LIST_PATH}" | tr -d '\r' | sed '/^$/d')
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
# Perform hadd for each run
#==============================

do_hadd() {
    local run="$1"

    local bucket_val=$(( (10#$run / 1000) * 1000 ))
    local group_val=$(( (10#$run / 100) * 100 ))
    
    local run_bucket=$(printf "%08d" "$bucket_val")
    local run_group=$(printf "%08d" "$group_val")

    local base_dir="${OUTPUT_DIR_BASE}/${run_bucket}"

    local group_candidates=($(find "$base_dir" -maxdepth 1 -type d -regex ".*/${run_group}\(_v[0-9]+\)?$" | sort))

    if (( ${#group_candidates[@]} == 0 )); then
        log WARN "No matching run_group directory found under $base_dir for group $run_group"
        return
    fi

    local selected_group_dir="${group_candidates[-1]}"
    local run_dir="${selected_group_dir}/${run}"

    if [[ ! -d "$run_dir" ]]; then
        log WARN "Run directory not found: $run_dir"
        return
    fi

    local files=($(printf "%s\n" "$run_dir"/*.root | sort -V))
    local n_files=${#files[@]}
    if (( n_files == 0 )); then
        log WARN "No ROOT files to merge for run $run"
        return
    fi


    local output_dir="${OUTPUT_DIR_BASE}/summary"
    mkdir -p "${output_dir}"
    local output_file="${output_dir}/RUN.${run}.${OUTPUT_SUFFIX}"

    log INFO "Merging ${n_files} files in $run_dir ..."
    hadd -f "$output_file" "${files[@]}" && \
        log INFO "Merged successfully -> $output_file" || \
        log ERROR "hadd failed for run $run"
}

#==============================
# Main
#==============================

main() {
    parse_args "$@"
    load_run_list
    filter_runs

    for run in "${RUN_LIST[@]}"; do
        do_hadd "$run"
    done

    log INFO "All hadd operations completed."
}

main "$@"