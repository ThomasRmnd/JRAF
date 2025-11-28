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
FILE_RANGE=100
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
  --file-range <num>             Max number of consecutive files per job (default: $FILE_RANGE)
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
    log INFO "Building grouped file ranges (max ${FILE_RANGE} per group)..."

    RANGES=()
    local range_start=""
    local prev_idx=""

    for i in "${!RTRAW_LIST[@]}"; do
        f="${RTRAW_LIST[$i]}"
        fname=${f##*/}

        if [[ $fname =~ \.[0-9]{14}\.([0-9]+)_ ]]; then
            num="${BASH_REMATCH[1]}"
            num=$((10#$num))  # strip leading zeros
        else
            log WARN "Could not extract file number from: $fname"
            continue
        fi

        if [[ -z "$range_start" ]]; then
            range_start=$i
            count=1
        else
            if (( num == prev_num + 1 && count < FILE_RANGE )); then
                ((count++))
            else
                RANGES+=("$range_start-$prev_idx")
                range_start=$i
                count=1
            fi
        fi
        prev_num=$num
        prev_idx=$i
    done

    RANGES+=("$range_start-$prev_idx")
    log INFO "Generated ${#RANGES[@]} job ranges"
}

#==============================
# Build extra args list
#==============================

prepare_extra_args() {
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
    for r in "${RANGES[@]}"; do
        start=${r%-*}
        end=${r#*-}

        log INFO "Submitting job for run ${RUN_NUMBER} range ${start}-${end}"

        sbatch \
            --job-name="agrpc_${RUN_NUMBER}_${start}_${end}_batch" \
            --output="/dev/null" \
            --error="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/err/agrpc_${RUN_NUMBER}_${start}_${end}.err" \
            --partition="htc" \
            --ntasks=1 \
            --cpus-per-task=1 \
            --mem="4G" \
            --time="0-08:00:00" \
            --mail-user="thomas.raymond@iphc.cnrs.fr" \
            --mail-type="FAIL" \
            job_worker.sh \
            "$RUN_NUMBER" "$start" "$end" --no-local-copy "${EXTRA_ARGS[@]}"
        # "/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/log/agrpc_${RUN_NUMBER}_${start}_${end}.log"
    done

    log INFO "All jobs submitted successfully"
}

#==============================
# Main
#==============================

main() {
    parse_args "$@"
    load_file_lists
    prepare_job_arrays
    prepare_extra_args
    submit_jobs
}

main "$@"