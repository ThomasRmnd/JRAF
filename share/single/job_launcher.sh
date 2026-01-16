#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO Job Submission Helper
#  Purpose: Automate job submissions for single processing
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

LIST_BASE_REPROD25C="/eos/juno/groups/DataQuality/P25A/Physics/goodrunlist_v3.6"
LIST_BASE_REPROD25D="/eos/juno/groups/DataQuality/ReProd25D/Physics/goodrunlist_v0.0"

TIME_WINDOW=("-2.0" "2.0")
LOG_LEVEL=3

#==============================
# Usage & Argument Parsing
#==============================

usage() {
    cat <<EOF
Usage: $(basename "$0") --site <str> --campaign <str> --run <number> [options]

Required:
  --site <str>                 Storage site selection {EOS|CNAF}
  --campaign <str>             Campaign selection {Normal|ReProd25A|ReProd25B|ReProd25C|ReProd25D}
  --run <num>                  Run number to process     

Optional:
  --range <num>                Number of files to process (default: all)
  --property-file <path>       Path to property file
  --time-window <min> <max>    Time window (default: ${TIME_WINDOW[*]})
  --log-level <num>            Logging level (default: $LOG_LEVEL)
  --help                       Show this help message and exit
EOF
}

parse_args() {
    if (( $# < 6 )); then
        echo "ERROR: Missing arguments" >&2
        usage >&2
        exit 1
    fi

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --site)          SITE="$2"; shift 2 ;;
            --campaign)      CAMPAIGN="$2"; shift 2 ;;
            --run)           RUN_NUMBER="$2"; shift 2 ;;
            --range)         RANGE="$2"; shift 2 ;;
            --property-file) PROPERTY_FILE="$2"; shift 2 ;;
            --time-window)   TIME_WINDOW=("$2" "$3"); shift 3 ;;
            --log-level)     LOG_LEVEL="$2"; shift 2 ;;
            --help|-h)       usage; exit 0 ;;
            *) log ERROR "Unknown argument: $1"; usage; exit 1 ;;
        esac
    done

    if [[ -z "${RUN_NUMBER:-}" ]]; then
        log ERROR "--run is required"
        usage
        exit 1
    fi

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

    if [[ -z "${CAMPAIGN:-}" ]]; then
        log ERROR "--campaign is required {Normal|ReProd25A|ReProd25B|ReProd25C|ReProd25D}"
        usage
        exit 1
    fi

    case "${CAMPAIGN}" in
        Normal)
            LIST_BASE="${LIST_BASE_REPROD25C}"
            ;;
        ReProd25A)
            LIST_BASE="${RUN_LIST_REPROD25C}"
            ;;
        ReProd25B)
            LIST_BASE="${LIST_BASE_REPROD25C}"
            ;;
        ReProd25C)
            LIST_BASE="${LIST_BASE_REPROD25C}"
            ;;
        ReProd25D)
            LIST_BASE="${LIST_BASE_REPROD25D}"
            ;;
        *)
            log ERROR "Invalid --campaign: ${CAMPAIGN}"
            usage
            exit 1
            ;;
    esac

    if [[ "${CLUSTER}" == "IHEP" && "${SITE}" == "CNAF" ]]; then
        log WARN "CNAF site was selected while running on IHEP cluster"
    fi
}

#==============================
# Load and validate file lists
#==============================

load_file_lists() {
    local rtraw_list_file="${LIST_BASE}/rtraw_list/run_${RUN_NUMBER}.txt"
    local esd_list_file="${LIST_BASE}/esd_list/run_${RUN_NUMBER}.txt"

    log INFO "Listing ROOT files from EOS..."
    mapfile -t RTRAW_LIST < <(xrdfs "${XRD_URL_EOS}" cat "${rtraw_list_file}")
    mapfile -t ESD_LIST   < <(xrdfs "${XRD_URL_EOS}" cat "${esd_list_file}")

    log INFO "Number of RTRAW files: ${#RTRAW_LIST[@]}"
    log INFO "Number of ESD   files: ${#ESD_LIST[@]}"
}

#==============================
# Prepare submission arguments
#==============================

prepare_job_arrays() {
    RANGE="${RANGE:-${#RTRAW_LIST[@]}}"

    if ! [[ "${RANGE}" =~ ^[0-9]+$ ]]; then
        log ERROR "file-offset and file-range must be non-negative integers"
        exit 1
    fi

    log INFO "Processing range: count=${RANGE}"

    RTRAW_LIST=("${RTRAW_LIST[@]:0:$RANGE}")
    ESD_LIST=("${ESD_LIST[@]:0:$RANGE}")

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

    if [[ "${CLUSTER}" == "CC-IN2P3" ]]; then
        PROPERTY_FILE="${PROPERTY_FILE:-/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/other/Reconstruction/RecMuon/CdWpTtChi2RecTool/CdWpTtChi2RecTool.Properties.json}"
    elif [[ "${CLUSTER}" == "IHEP" ]]; then
        PROPERTY_FILE="${PROPERTY_FILE:-/junofs/users/traymond/reconstruction/esd/properties/RUN.${RUN_NUMBER}.Properties.json}"
    else
        log ERROR "No property file suitable in ${CLUSTER} cluster"
        exit 1
    fi

    EXTRA_ARGS=(
        "--property-file" "${PROPERTY_FILE}"
        "--time-window" "${TIME_WINDOW[0]}" "${TIME_WINDOW[1]}"
        "--log-level" "${LOG_LEVEL}"
    )
}

#==============================
# Submit jobs
#==============================

select_slurm_time_limit() {
    case "${SITE}" in
        EOS)
            SBATCH_TIME="0-00:30:00"
            ;;
        CNAF)
            SBATCH_TIME="0-00:20:00"
            ;;
        *)
            log ERROR "Unknown site for SLURM time selection: ${SITE}"
            exit 1
            ;;
    esac
}

submit_jobs() {
    log INFO "Submitting ${JOB_COUNT} jobs (cluster=${CLUSTER}, site=${SITE}, campaign=${CAMPAIGN})..."
    if [[ "${CLUSTER}" == "CC-IN2P3" ]]; then
        #------------------------------
        #  SLURM submission (IN2P3)
        #------------------------------
        select_slurm_time_limit
        
        sbatch \
            --job-name="agrpc_${RUN_NUMBER}_batch" \
            --output="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/log/agrpc_${RUN_NUMBER}_%a.log" \
            --error="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/err/agrpc_${RUN_NUMBER}_%a.err" \
            --array="0-$((JOB_COUNT - 1))" \
            --partition="htc" \
            --ntasks=1 \
            --cpus-per-task=1 \
            --mem="3G" \
            --time="${SBATCH_TIME}" \
            --mail-user="thomas.raymond@iphc.cnrs.fr" \
            --mail-type="FAIL" \
            job_worker.sh \
            "${RUN_NUMBER}" "${SITE}" "${CAMPAIGN}" "${EXTRA_ARGS[@]}"

    elif [[ "${CLUSTER}" == "IHEP" ]]; then
        #------------------------------
        #  hep_sub submission (IHEP)
        #------------------------------

        hep_sub job_worker.sh \
            -argu "%{ProcId} ${RUN_NUMBER} ${SITE} ${CAMPAIGN} ${EXTRA_ARGS[*]}" \
            -n "$JOB_COUNT" \
            -cpu 1 \
            -m 4096 \
            -wt short \
            -o "/scratchfs/juno/traymond/agrpc_${RUN_NUMBER}_%{ProcId}.log" \
            -e "/scratchfs/juno/traymond/agrpc_${RUN_NUMBER}_%{ProcId}.err" \
            -name "agrpc_${RUN_NUMBER}_batch"

    else
        log ERROR "Unknown cluster: ${CLUSTER}. Cannot submit jobs."
        exit 1
    fi
    log INFO "Jobs submitted successfully"
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