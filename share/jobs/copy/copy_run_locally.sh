#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO Run Downloader
#  Purpose: Copy all ReProd files for a specific run from XRD to local storage
#--------------------------------------------------------------------------------------------------

set -euo pipefail
IFS=$'\n\t'

#==============================
#  Utility logging function
#==============================

RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
BLUE='\033[0;34m'
NC='\033[0m'

log() {
    local level="$1"; shift
    local msg="$*"
    local timestamp
    timestamp="$(date '+%Y-%m-%d %H:%M:%S')"

    local level_num=0 color="$NC"
    case "$level" in
        ERROR) level_num=1; color="$RED" ;;
        WARN)  level_num=2; color="$YELLOW" ;;
        INFO)  level_num=3; color="$GREEN" ;;
        DEBUG) level_num=4; color="$CYAN" ;;
        ALL)   level_num=5; color="$BLUE" ;;
        *)     level_num=3 ;;
    esac

    local prefix="${color}[$timestamp][$level]${NC}"

    case "$level" in
        DEBUG|INFO) echo -e "${prefix} $msg" >&1 ;;
        WARN|ERROR) echo -e "${prefix} $msg" >&2 ;;
        ALL)
            echo -e "${prefix} $msg" >&1
            echo -e "${prefix} $msg" >&2
            ;;
        *) echo -e "${prefix} $msg" >&1 ;;
    esac
}

#==============================
# Utility functions & Cluster Detection
#==============================

HOSTNAME=$(hostname -f 2>/dev/null || hostname)
if [[ "${HOSTNAME}" =~ ^cca[0-9]+\.in2p3\.fr$ ]]; then
    CLUSTER="CC-IN2P3"
    export X509_USER_PROXY=/sps/juno/jdeandre/rtraw_ThomasRaymond/.cert_traymond_juno_user
elif [[ "${HOSTNAME}" =~ ^lxlogin[0-9]+\.ihep\.ac\.cn$ ]]; then
    CLUSTER="IHEP"
else
    echo "ERROR: Unknown cluster: ${HOSTNAME}" >&2
    exit 1
fi

log INFO "Cluster detected: ${CLUSTER}"

#==============================
# Configuration
#==============================

RUN_LIST_REPROD25C="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/other/GoodList/ReProd25C/physics_good.txt"
RUN_LIST_REPROD25D="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/other/GoodList/ReProd25D/physics_good.txt"

XRD_URL_EOS="root://junoeos01.ihep.ac.cn/"
XRD_URL_CNAF="root://xrootd-archive.cr.cnaf.infn.it:1095/"

XRD_BASEPATH_EOS="/eos"
XRD_BASEPATH_CNAF="/production/storm/dirac"

LOCAL_CCIN2P3_BASEPATH="/sps/juno/jdeandre/rtraw_ThomasRaymond"
LOCAL_IHEP_BASEPATH="/junofs/users/traymond"

#==============================
# Usage & Argument Parsing
#==============================

usage() {
    cat <<EOF
Usage: $(basename "$0") --site <str> --campaign <str> --run <int> [options]

Required:
  --site <str>          Remote site {EOS|CNAF}
  --campaign <str>      Campaign {Normal|ReProd25A|ReProd25B|ReProd25C|ReProd25D}
  --run <int>           Run number to download

Optional:
  --stream <str>        Data stream (default: global_trigger)
  --force               Overwrite existing local files
EOF
}

parse_args() {
    STREAM="global_trigger"
    FORCE=0

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --site)     SITE="$2"; shift 2 ;;
            --campaign) CAMPAIGN="$2"; shift 2 ;;
            --run)      RUN_NUMBER="$2"; shift 2 ;;
            --stream)   STREAM="$2"; shift 2 ;;
            --help|-h)  usage; exit 0 ;;
            *) log ERROR "Unknown argument: $1"; usage; exit 1 ;;
        esac
    done

    [[ -z "${SITE:-}" ]] && { log ERROR "--site is required"; exit 1; }
    [[ -z "${CAMPAIGN:-}" ]] && { log ERROR "--campaign is required"; exit 1; }
    [[ -z "${RUN_NUMBER:-}" ]] && { log ERROR "--run is required"; exit 1; }

    case "${CLUSTER}" in
        IHEP)
            LOCAL_BASEPATH="${LOCAL_IHEP_BASEPATH}"
            ;;
        CC-IN2P3)
            LOCAL_BASEPATH="${LOCAL_CCIN2P3_BASEPATH}"
            ;;
        *)
            log ERROR "Unknown cluster: ${CLUSTER}"
            exit 1
            ;;
    esac

    case "${SITE}" in
        EOS)
            XRD_URL="${XRD_URL_EOS}"
            XRD_BASEPATH="${XRD_BASEPATH_EOS}"
            ;;
        CNAF)
            XRD_URL="${XRD_URL_CNAF}"
            XRD_BASEPATH="${XRD_BASEPATH_CNAF}"
            PROXY_PATH="/sps/juno/jdeandre/rtraw_ThomasRaymond/.cert_traymond_juno_user"
            if [[ ! -f "${PROXY_PATH}" ]]; then
                log ERROR "X.509 proxy does not exist: ${PROXY_PATH}"
                exit 1
            fi
            if [[ ! -r "${PROXY_PATH}" ]]; then
                log ERROR "X.509 proxy not readable: ${PROXY_PATH}"
                exit 1
            fi
            export X509_USER_PROXY="${PROXY_PATH}"
            ;;
        *)
            log ERROR "Invalid site argument: ${SITE} (expected {EOS|CNAF})"
            exit 1
            ;;
    esac

    case "${CAMPAIGN}" in
        Normal)
            LIST_BASE="${RUN_LIST_REPROD25C%/*}"
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            ;;
        ReProd25A)
            LIST_BASE="${RUN_LIST_REPROD25C%/*}"
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            ;;
        ReProd25B)
            LIST_BASE="${RUN_LIST_REPROD25C%/*}"
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            ;;
        ReProd25C)
            LIST_BASE="${RUN_LIST_REPROD25C%/*}"
            RUN_LIST_PATH="${RUN_LIST_REPROD25C}"
            ;;
        ReProd25D)
            LIST_BASE="${RUN_LIST_REPROD25D%/*}"
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
# Path Resolution
#==============================

resolve_remote_directory() {
    local bucket_val=$(( (10#$RUN_NUMBER / 1000) * 1000 ))
    local group_val=$(( (10#$RUN_NUMBER / 100) * 100 ))
    local run_bucket=$(printf "%08d" "$bucket_val")
    local run_group=$(printf "%08d" "$group_val")

    local base_search_dir="${XRD_BASEPATH}/juno/juno-reprod/${CAMPAIGN}/${STREAM}/${run_bucket}"
    
    log INFO "Searching for latest group version in ${base_search_dir}..."

    local remote_group
    remote_group=$(xrdfs "${XRD_URL}" ls "${base_search_dir}" 2>/dev/null | grep -E "/${run_group}(_v[0-9]+)?/?$" | sort | tail -n 1)

    if [[ -z "${remote_group}" ]]; then
        log ERROR "Could not find group ${run_group} for run ${RUN_NUMBER} on ${SITE}"
        exit 1
    fi

    REMOTE_DIR="${remote_group}/${RUN_NUMBER}"

    local group_name=$(basename "${remote_group}")
    LOCAL_DIR="${LOCAL_BASEPATH}/juno/juno-reprod/${CAMPAIGN}/${STREAM}/${run_bucket}/${group_name}/${RUN_NUMBER}"
}

#==============================
# Execution
#==============================

download_run() {
    log INFO "Remote Directory: ${XRD_URL}${REMOTE_DIR}"
    log INFO "Local Directory:  ${LOCAL_DIR}"

    mkdir -p "${LOCAL_DIR}"

    mapfile -t FILE_LIST < <(xrdfs "${XRD_URL}" ls "${REMOTE_DIR}" | grep "RUN\.${RUN_NUMBER}.*\.esd$")

    if (( ${#FILE_LIST[@]} == 0 )); then
        log WARN "No ROOT files found for run ${RUN_NUMBER}."
        return
    fi

    log INFO "Found ${#FILE_LIST[@]} files. Starting download..."

    for remote_file in "${FILE_LIST[@]}"; do
        local filename=$(basename "${remote_file}")
        local dest_file="${LOCAL_DIR}/${filename}"

        if [[ -f "${dest_file}" && $FORCE -eq 0 ]]; then
            log DEBUG "File ${filename} already exists locally, skipping."
            continue
        fi

        log INFO "Copying ${filename}..."
        if xrdcp -f "${XRD_URL}${remote_file}" "${dest_file}"; then
            log INFO "Success: ${filename}"
        else
            log ERROR "Failed to copy ${filename}"
        fi
    done
}

#==============================
# Main
#==============================

main() {
    parse_args "$@"
    resolve_remote_directory
    download_run
    log INFO "Operation completed for run ${RUN_NUMBER}."
}

main "$@"