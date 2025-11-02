#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO Job Submission Helper
#  Purpose: Automate hep_sub job submissions for single ESD–RTRAW processing
#--------------------------------------------------------------------------------------------------

set -euo pipefail
IFS=$'\n\t'

#==============================
# Utility functions
#==============================

source /junofs/users/traymond/bash/logging.sh

#==============================
# Configuration defaults
#==============================

EOS_BASE="root://junoeos01.ihep.ac.cn/"
LIST_BASE="/eos/juno/groups/DataQuality/P25A/Physics/goodrunlist_v3.4"
TIME_WINDOW=("-2.0" "2.0")
LOG_LEVEL=3

#==============================
# Usage & Argument Parsing
#==============================

usage() {
    cat <<EOF
Usage: $(basename "$0") --run-number <number> [options]

Required:
  --run-number <number>        Run number to process

Optional:
  --list-base <path>           Path to list base (default: $LIST_BASE)
  --file-offset <num>          Starting index in the file list (default: 0)
  --file-range <num>           Number of files to process (default: all)
  --property-file <path>       Path to property file
  --time-window <min> <max>    Time window (default: ${TIME_WINDOW[*]})
  --log-level <num>            Logging level (default: $LOG_LEVEL)
  --help                       Show this help message and exit
EOF
}

parse_args() {
    if [[ $# -eq 0 ]]; then
        usage
        exit 1
    fi

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --run-number)    RUN_NUMBER="$2"; shift 2 ;;
            --list-base)     LIST_BASE="$2"; shift 2 ;;
            --file-offset)   FILE_OFFSET="$2"; shift 2 ;;
            --file-range)    FILE_RANGE="$2"; shift 2 ;;
            --property-file) PROPERTY_FILE="$2"; shift 2 ;;
            --time-window)   TIME_WINDOW=("$2" "$3"); shift 3 ;;
            --log-level)     LOG_LEVEL="$2"; shift 2 ;;
            --help|-h)       usage; exit 0 ;;
            *) log ERROR "Unknown argument: $1"; usage; exit 1 ;;
        esac
    done

    if [[ -z "${RUN_NUMBER:-}" ]]; then
        log ERROR "--run-number is required"
        usage
        exit 1
    fi
}

#==============================
# Load and validate file lists
#==============================

load_file_lists() {
    local rtraw_list_file="${LIST_BASE}/rtraw_list/run_${RUN_NUMBER}.txt"
    local esd_list_file="${LIST_BASE}/esd_list/run_${RUN_NUMBER}.txt"

    log INFO "Listing ROOT files from EOS..."
    mapfile -t RTRAW_LIST < <(xrdfs "$EOS_BASE" cat "$rtraw_list_file")
    mapfile -t ESD_LIST   < <(xrdfs "$EOS_BASE" cat "$esd_list_file")

    log INFO "Number of RTRAW files: ${#RTRAW_LIST[@]}"
    log INFO "Number of ESD   files: ${#ESD_LIST[@]}"
}


#==============================
# Prepare submission arguments
#==============================

prepare_job_arrays() {
    FILE_OFFSET="${FILE_OFFSET:-0}"
    FILE_RANGE="${FILE_RANGE:-$(( ${#RTRAW_LIST[@]} - FILE_OFFSET ))}"

    if ! [[ "$FILE_OFFSET" =~ ^[0-9]+$ && "$FILE_RANGE" =~ ^[0-9]+$ ]]; then
        log ERROR "file-offset and file-range must be non-negative integers"
        exit 1
    fi

    log INFO "Processing range: offset=$FILE_OFFSET, count=$FILE_RANGE"

    RTRAW_LIST=("${RTRAW_LIST[@]:$FILE_OFFSET:$FILE_RANGE}")
    ESD_LIST=("${ESD_LIST[@]:$FILE_OFFSET:$FILE_RANGE}")

    JOB_COUNT_RTRAW=${#RTRAW_LIST[@]}
    JOB_COUNT_ESD=${#ESD_LIST[@]}

    if (( JOB_COUNT_RTRAW != JOB_COUNT_ESD )); then
        log ERROR "Mismatch between RTRAW ($JOB_COUNT_RTRAW) and ESD ($JOB_COUNT_ESD) files"
        exit 1
    fi

    JOB_COUNT=${JOB_COUNT_RTRAW}

    if (( JOB_COUNT == 0 )); then
        log ERROR "No ROOT files found for run ${RUN_NUMBER}"
        exit 1
    fi

    PROPERTY_FILE="${PROPERTY_FILE:-/junofs/users/traymond/reconstruction/esd/properties/RUN.${RUN_NUMBER}.Properties.json}"

    EXTRA_ARGS=(
        "--property-file" "$PROPERTY_FILE"
        "--time-window" "${TIME_WINDOW[0]}" "${TIME_WINDOW[1]}"
        "--log-level" "$LOG_LEVEL"
    )
}

#==============================
# Submit jobs via hep_sub
#==============================

submit_jobs() {
    log INFO "Submitting ${JOB_COUNT} jobs via hep_sub..."

    hep_sub job_worker.sh \
        -argu "%{ProcId} ${RUN_NUMBER} ${LIST_BASE} ${EXTRA_ARGS[*]}" \
        -n "$JOB_COUNT" \
        -cpu 1 \
        -m 4096 \
        -wt short \
        -o "/scratchfs/juno/traymond/agrpc_${RUN_NUMBER}_%{ProcId}.log" \
        -e "/scratchfs/juno/traymond/agrpc_${RUN_NUMBER}_%{ProcId}.err" \
        -name "agrpc_${RUN_NUMBER}_batch"

    log INFO "All jobs submitted successfully"
}

#==============================
# Main
#==============================

main() {
    parse_args "$@"
    load_file_lists
    prepare_job_arrays
    submit_jobs
}

main "$@"