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
  --lower <num>         Starting run number (inclusive)
  --upper <num>         Ending run number (inclusive)
  --help           Show this help message and exit

Examples:
  $(basename "$0") --lower 10500 --upper 10550
EOF
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --lower)   LOWER_BOUND="$2"; shift 2 ;;
            --upper)   UPPER_BOUND="$2"; shift 2 ;;
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
    mapfile -t RUN_LIST < <(xrdfs "${EOS_BASE}" cat "${RUN_LIST_PATH}" | tr -d '\r' | sed '/^$/d')
}

#==============================
# Range Filtering
#==============================

filter_runs() {
    local filtered=()

    for run in "${RUN_LIST[@]}"; do
        (( run < 0 )) && continue
        if [[ -n "${LOWER_BOUND}" && "$run" -lt "${LOWER_BOUND}" ]]; then
            continue
        fi
        if [[ -n "${UPPER_BOUND}" && "$run" -gt "${UPPER_BOUND}" ]]; then
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
        sbatch \
            --job-name="comp_${run}_batch" \
            --output="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/log/comp_${run}.log" \
            --error="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/err/comp_${run}.err" \
            --partition="htc" \
            --ntasks=1 \
            --cpus-per-task=1 \
            --mem="8G" \
            --time="0-01:00:00" \
            --mail-user="thomas.raymond@iphc.cnrs.fr" \
            --mail-type="FAIL" \
            job_worker.sh \
            "${run}"
        # "/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/log/comp_${run}_%a.log"
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