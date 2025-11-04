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

source /pbs/home/t/traymond/share/bash/logging.sh

#==============================
# Configuration defaults
#==============================

EOS_BASE="root://junoeos01.ihep.ac.cn/"
BASE_OUTPUT_DIR="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/ibd"
RUN_LIST_PATH="/eos/juno/groups/DataQuality/P25A/Physics/goodrunlist_v3.4/Physics_good_run_list.txt"

LOWER_BOUND=""
UPPER_BOUND=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --lower <num>      Starting run number (inclusive)
  --upper <num>      Ending run number (inclusive)
  --list  <path>     Path to custom run list (default: $RUN_LIST_PATH)
  --help             Show this help message and exit

Examples:
  $(basename "$0") --lower 9000 --upper 9050
  $(basename "$0") --list /path/to/custom_run_list.txt
EOF
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --lower) LOWER_BOUND="$2"; shift 2 ;;
            --upper) UPPER_BOUND="$2"; shift 2 ;;
            --list)  RUN_LIST_PATH="$2"; shift 2 ;;
            --help|-h) usage; exit 0 ;;
            *) log ERROR "Unknown argument: $1"; usage; exit 1 ;;
        esac
    done
}

#==============================
# Load and filter run list
#==============================

load_run_list() {
    log INFO "Fetching run list..."
    mapfile -t RUN_LIST < <(xrdfs "$EOS_BASE" cat "$RUN_LIST_PATH" | tr -d '\r' | sed '/^$/d')
}

filter_runs() {
    local filtered=()
    for run in "${RUN_LIST[@]}"; do
        (( run < 0 )) && continue
        if [[ -n "$LOWER_BOUND" && "$run" -lt "$LOWER_BOUND" ]]; then
            continue
        fi
        if [[ -n "$UPPER_BOUND" && "$run" -gt "$UPPER_BOUND" ]]; then
            continue
        fi
        filtered+=("$run")
    done
    RUN_LIST=("${filtered[@]}")

    if (( ${#RUN_LIST[@]} == 0 )); then
        log WARN "No runs matched range."
        exit 0
    fi

    log INFO "Runs to hadd: ${RUN_LIST[*]}"
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

    local base_dir="${BASE_OUTPUT_DIR}/${run_bucket}"

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

    local output_dir="${BASE_OUTPUT_DIR}/summary"
    mkdir -p "$output_dir"
    local output_file="${output_dir}/RUN.${run}.output.reprod.cca.root"

    log INFO "Merging files in $run_dir ..."
    local files=("$run_dir"/*.root)

    if (( ${#files[@]} == 0 )); then
        log WARN "No ROOT files to merge for run $run"
        return
    fi

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
