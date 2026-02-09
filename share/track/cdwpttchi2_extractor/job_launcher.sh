#!/bin/bash
#--------------------------------------------------------------------------------------------------
#  JUNO Multi-job Launcher
#  Purpose: Automate job multi-submissions
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
RUN_LIST_REPROD25D="/eos/juno/groups/DataQuality/ReProd25D/Physics/goodrunlist_v0.0/physics_good_v0.0.txt"

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
# Launch Jobs
#==============================

launch_jobs() {
    for run in "${RUN_LIST[@]}"; do
        log INFO ">>> Launching job for run ${run}"
        sbatch \
            --job-name="cdwptt_${run}_batch" \
            --output="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/log/cdwptt_${run}.log" \
            --error="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/err/cdwptt_${run}.err" \
            --partition="htc" \
            --ntasks=1 \
            --cpus-per-task=1 \
            --mem="2G" \
            --time="0-00:15:00" \
            --mail-user="thomas.raymond@iphc.cnrs.fr" \
            --mail-type="FAIL" \
            job_worker.sh \
            "${CAMPAIGN}" "${run}" "${LIST_BASE}"
        # "/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/log/cdwptt_${run}_%a.log"
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