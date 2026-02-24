#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO Job Submission Helper
#  Purpose: Automate job submissions for multi processing
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

TIME_WINDOW=("-2.0" "2.0")
LOG_LEVEL=3

#==============================
# Usage & Argument Parsing
#==============================

usage() {
    cat <<EOF
Usage: $(basename "$0") --iste <str> --campaign <str> --run <int> --list-base <str> --range <int> [options]

Required:
  --site <str>                  Storage site selection {EOS|CNAF}
  --campaign <str>              Campaign selection {Normal|ReProd25A|ReProd25B|ReProd25C|ReProd25D}
  --run <int>                   Run number to process
  --list-base <str>             Basepath for the file list
  --range <int>                 Number of files to process per job

Optional:
  --property-file <path>        Path to property file
  --time-window <float> <float> Time window (default: ${TIME_WINDOW[*]})
  --log-level <int>             Logging level (default: $LOG_LEVEL)
  --help                        Show this help message and exit
EOF
}

parse_args() {
    if [[ $# -eq 0 ]]; then
        usage
        exit 1
    fi

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --site)          SITE="$2"; shift 2 ;;
            --campaign)      CAMPAIGN="$2"; shift 2 ;;
            --run)           RUN_NUMBER="$2"; shift 2 ;;
            --list-base)     LIST_BASE="$2"; shift 2 ;;
            --range)         RANGE="$2"; shift 2 ;;
            --property-file) PROPERTY_FILE="$2"; shift 2 ;;
            --time-window)   TIME_WINDOW=("$2" "$3"); shift 3 ;;
            --log-level)     LOG_LEVEL="$2"; shift 2 ;;
            --help|-h)       usage; exit 0 ;;
            *) log ERROR "Unknown argument: $1"; usage; exit 1 ;;
        esac
    done

    if [[ -z "${SITE:-}" ]]; then
        log ERROR "--site is required {EOS|CNAF}"
        usage
        exit 1
    fi

    case "${SITE}" in
        EOS|CNAF) ;;
        *) log ERROR "Invalid --site: ${SITE} (expected {EOS|CNAF})"
           exit 1 ;;
    esac

    if [[ "${CLUSTER}" == "IHEP" && "${SITE}" == "CNAF" ]]; then
        log WARN "CNAF site was selected while running on IHEP cluster"
    fi

    if [[ -z "${CAMPAIGN:-}" ]]; then
        log ERROR "--campaign is required {Normal|ReProd25A|ReProd25B|ReProd25C|ReProd25D}"
        usage
        exit 1
    fi

    case "${CAMPAIGN}" in
        Normal|ReProd25A|ReProd25B|ReProd25C|ReProd25D) ;;
        *) log ERROR "Invalid --site: ${CAMPAIGN} (expected {Normal|ReProd25A|ReProd25B|ReProd25C|ReProd25D})"
           exit 1 ;;
    esac

    if [[ -z "${RUN_NUMBER:-}" ]]; then
        log ERROR "--run is required"
        usage
        exit 1
    fi

    if [[ -z "${LIST_BASE:-}" ]]; then
        log ERROR "--list-base is required"
        usage
        exit 1
    fi

    if [[ -z "${RANGE:-}" ]]; then
        log ERROR "--range is required"
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
    mapfile -t RTRAW_LIST < <(cat "${rtraw_list_file}")
    mapfile -t ESD_LIST   < <(cat "${esd_list_file}")

    log INFO "Number of RTRAW files: ${#RTRAW_LIST[@]}"
    log INFO "Number of ESD   files: ${#ESD_LIST[@]}"
}

#==============================
# Prepare submission arguments
#==============================

prepare_job_arrays() {
    log INFO "Building grouped file ranges (max ${RANGE} per group)..."

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
            if (( num == prev_num + 1 && count < RANGE )); then
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
    PROPERTY_FILE="${PROPERTY_FILE:-/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/other/Reconstruction/RecMuon/CdWpTtChi2RecTool/CdWpTtChi2RecTool.Properties.json}"

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
            --output="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/log/agrpc_${RUN_NUMBER}_${start}_${end}.log" \
            --error="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/err/agrpc_${RUN_NUMBER}_${start}_${end}.err" \
            --partition="htc" \
            --ntasks=1 \
            --cpus-per-task=1 \
            --mem="4G" \
            --time="0-8:00:00" \
            --mail-user="thomas.raymond@iphc.cnrs.fr" \
            --mail-type="FAIL" \
            job_worker.sh \
            "${SITE}" "${CAMPAIGN}" "${RUN_NUMBER}" "${LIST_BASE}" "${start}" "${end}" "${EXTRA_ARGS[@]}"
        # --local --no-local-copy --skip-if-exist
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