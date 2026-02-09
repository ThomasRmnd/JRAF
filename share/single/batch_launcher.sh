#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO Job Multi-Submission Helper
#  Purpose: Automate job multi-submissions for single processing
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

LOWER_BOUND=""
UPPER_BOUND=""

usage() {
    cat <<EOF
Usage: $(basename "$0") --site <str> --campaign <str> [options]

Required:
  --site <str>          Storage site selection {EOS|CNAF}
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
            --site)    SITE="$2"; shift 2 ;;
            --campaign) CAMPAIGN="$2"; shift 2 ;;
            --lower)   LOWER_BOUND="$2"; shift 2 ;;
            --upper)   UPPER_BOUND="$2"; shift 2 ;;
            --help|-h) usage; exit 0 ;;
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
        *) 
            log ERROR "Invalid --site: ${SITE}"
            usage
            exit 1 ;;
    esac

    if [[ "${CLUSTER}" == "IHEP" && "${SITE}" == "CNAF" ]]; then
        log WARN "CNAF site was selected while running on IHEP cluster"
    fi

    if [[ -z "${CAMPAIGN:-}" ]]; then
        log ERROR "--campaign is required"
        usage
        exit 1
    fi

    case "${CAMPAIGN}" in
        Normal)
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            ;;
        ReProd25A)
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            ;;
        ReProd25B)
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            ;;
        ReProd25C)
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            ;;
        ReProd25D)
            RUN_LIST_PATH="${RUN_LIST_REPROD25D}"
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
# Launch Jobs
#==============================

launch_jobs() {
    for run in "${RUN_LIST[@]}"; do
        log INFO ">>> Launching job for run ${run}"
        local cmd=(sh job_launcher.sh --site ${SITE} --campaign ${CAMPAIGN} --run ${run})

        if "${cmd[@]}"; then
            log INFO "Run ${run} submitted successfully"
        else
            log ERROR "Submission failed for run ${run}"
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