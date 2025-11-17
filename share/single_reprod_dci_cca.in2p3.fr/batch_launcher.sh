#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO Job Multi-Submission Helper
#  Purpose: Automate hep_sub job multi-submissions for ReProd processing
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
RUN_LIST_PATH="/eos/juno/groups/DataQuality/P25A/Physics/goodrunlist_v3.6/Physics_good_run_list.txt"

LOWER_BOUND=""
UPPER_BOUND=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  -r <num>         Starting run number (inclusive)
  -R <num>         Ending run number (inclusive)
  -l <path>        Path to custom run list (default: $RUN_LIST_PATH)
  --help           Show this help message and exit

Examples:
  $(basename "$0") --r 9000 --R 9050
  $(basename "$0") --list /path/to/custom_run_list.txt
EOF
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --lower)   LOWER_BOUND="$2"; shift 2 ;;
            --upper)   UPPER_BOUND="$2"; shift 2 ;;
            --list)    RUN_LIST_PATH="$2"; shift 2 ;;
            --help|-h) usage; exit 0 ;;
            *) log ERROR "Unknown argument: $1"; usage; exit 1 ;;
        esac
    done
}

#==============================
# Fetch Run List
#==============================

load_run_list() {
    log INFO "Fetching run list from EOS..."
    mapfile -t RUN_LIST < <(xrdfs "$EOS_BASE" cat "$RUN_LIST_PATH" | tr -d '\r' | sed '/^$/d')
}

#==============================
# Range Filtering
#==============================

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
        log WARN "No runs matched the provided range"
        exit 0
    fi

    log INFO "Runs selected: ${RUN_LIST[*]}"
}

#==============================
# Launch Jobs
#==============================

launch_jobs() {
    for run in "${RUN_LIST[@]}"; do
        log INFO ">>> Launching job for run ${run}"
        local cmd=(bash job_launcher.sh --run-number "$run")

        if "${cmd[@]}"; then
            log INFO "Run $run submitted successfully"
        else
            log ERROR "Submission failed for run $run"
        fi
    done
}

#==============================
# Main
#==============================

main() {
    parse_args "$@"
    load_run_list
    filter_runs
    launch_jobs
    log INFO "All matching runs processed successfully"
}

main "$@"