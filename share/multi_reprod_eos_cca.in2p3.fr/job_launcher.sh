#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO Job Submission Helper
#  Purpose: Automate hep_sub job submissions for multi ReProd processing
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
LIST_BASE="/eos/juno/groups/DataQuality/P25A/Physics/goodrunlist_v3.6"
TIME_WINDOW=("-2.0" "2.0")
LOG_LEVEL=3
FILE_RANGE=100

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

    log INFO "Listing ROOT files from EOS..."
    mapfile -t RTRAW_LIST   < <(xrdfs "$EOS_BASE" cat "$rtraw_list_file")

    log INFO "Number of RTRAW files: ${#RTRAW_LIST[@]}"
}


#==============================
# Prepare submission arguments
#==============================

prepare_job_arrays() {
    if ! [[ "$FILE_RANGE" =~ ^[0-9]+$ ]]; then
        log ERROR "file-range must be non-negative integers"
        exit 1
    fi

    log INFO "Processing range: $FILE_RANGE"

    RANGES=()
    range_start=""
    prev_num=""
    count=0

    for f in "${RTRAW_LIST[@]}"; do
        fname=${f##*/}

        if [[ $fname =~ \.[0-9]{14}\.([0-9]+)_ ]]; then
            file_number=${BASH_REMATCH[1]}
            file_number=$((10#$file_number))  # strip leading zeros
        else
            log WARN "Warning: could not extract file number from $fname"
            continue
        fi

        if [[ -z "$range_start" ]]; then
            range_start=$file_number
            count=1
        else
            if (( file_number == prev_num + 1 && count < $FILE_RANGE )); then
                ((count++))
            else
                RANGES+=("$range_start-$prev_num")
                range_start=$file_number
                count=1
            fi
        fi
        prev_num=$file_number
    done

    RANGES+=("$range_start-$prev_num")
    JOB_COUNT=${#RANGES[@]}

    if (( JOB_COUNT == 0 )); then
        log ERROR "No valid contiguous ranges found for run ${RUN_NUMBER}"
        exit 1
    fi

    log INFO "Identified ${JOB_COUNT} contiguous ranges"

    PROPERTY_FILE="${PROPERTY_FILE:-/sps/juno/jdeandre/rtraw_ThomasRaymond/esd/properties/RUN.${RUN_NUMBER}.Properties.json}"

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

    for r in "${RANGES[@]}"; do
        START=${r%-*}
        END=${r#*-}
        log INFO "Submitting job for run ${RUN_NUMBER} range ${START}-${END}"
        hep_sub job_worker.sh \
            -argu "${START} ${END} ${RUN_NUMBER} ${LIST_BASE} ${EXTRA_ARGS[*]}" \
            -cpu 1 \
            -m 4096 \
            -wt short \
            -o "/dev/null" \
            -e "/scratchfs/juno/traymond/agrpc_${RUN_NUMBER}_${START}_${END}.err" \
            -name "agrpc_${RUN_NUMBER}_${START}_${END}_batch"
        # -o "/scratchfs/juno/traymond/agrpc_${RUN_NUMBER}_${START}_${END}.log" \
    done

    log INFO "All jobs submitted successfully"
}

main() {
    parse_args "$@"
    load_file_lists
    prepare_job_arrays
    submit_jobs
}

main "$@"