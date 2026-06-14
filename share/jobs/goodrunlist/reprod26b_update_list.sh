#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO ReProd26B Good Run List Updater
#  Purpose: Update the Good Run List of the ReProd25B with the given list version
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
# Utility functions
#==============================

HOSTNAME=$(hostname -f 2>/dev/null || hostname)
if [[ "${HOSTNAME}" =~ ^cc.*\.in2p3\.fr$ ]]; then
    # Detect CC-IN2P3 cluster
    CLUSTER="CC-IN2P3"
elif [[ "${HOSTNAME}" =~ ^lxlogin[0-9]+\.ihep\.ac\.cn$ ]]; then
    # Detect IHEP cluster
    CLUSTER="IHEP"
else
    echo "ERROR: Unknown cluster. Hostname: ${HOSTNAME}" >&2
    echo "Expected CC-IN2P3 (cca###) or IHEP (lxlogin###.ihep.ac.cn)" >&2
    exit 1
fi

log INFO "Cluster detected: ${CLUSTER}"

#==============================
# Safety helpers
#==============================

require_dir() {
    [[ -d "$1" ]] || {
        log ERROR "Missing directory: $1"
        exit 1
    }
}

safe_rm_file() {
    local path="$1"
    if [[ -e "${path}" ]]; then
        log WARN "Removing: ${path}"
        rm -- "${path}"
    fi
}

safe_rm_dir() {
    local path="$1"
    if [[ -d "${path}" ]]; then
        log WARN "Removing: ${path}"
        rm -r -- "${path}"
    fi
}

#==============================
# Configuration defaults
#==============================

GOODRUNLIST_PATH="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/other/GoodList"

ORIG_BASE="${GOODRUNLIST_PATH}/ReProd26B_Original/Physics"
DEST_BASE="${GOODRUNLIST_PATH}/ReProd26B"

PHASES=(phase1 phase2 phase3 phase4 phase4b)

#==============================
# Usage & Argument Parsing
#==============================

VERSION=""

usage() {
    cat <<EOF
Usage: $(basename "$0") --version <str>

Required:
  --version             <str>           Version of the Good Run List (example: v5.0.5)

Optional:
  --help                                Show this help message and exit
EOF
}

parse_args() {
    if [[ $# -eq 0 ]]; then
        usage
        exit 1
    fi

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --version)       VERSION="$2"; shift 2 ;;
            --help|-h)       usage; exit 0 ;;
            *) log ERROR "Unknown argument: $1"; usage; exit 1 ;;
        esac
    done

    if [[ -z "${VERSION:-}" ]]; then
        log ERROR "--version is required"
        usage
        exit 1
    fi

    if [[ ! "${VERSION}" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
        log ERROR "Invalid version format: '${VERSION}' (expected e.g. v5.0.5)"
        exit 1
    fi

    log INFO "Version: ${VERSION}"
}


#==============================
# Clean
#==============================

clean_dir() {
    require_dir "${DEST_BASE}"

    safe_rm_file "${DEST_BASE}/physics_good.txt"
    for ph in "${PHASES[@]}"; do
        safe_rm_file "${DEST_BASE}/physics_good_${ph}.txt"
    done

    safe_rm_dir "${DEST_BASE}/rtraw_list"
    safe_rm_dir "${DEST_BASE}/esd_list"

    mkdir -p "${DEST_BASE}/rtraw_list"
    mkdir -p "${DEST_BASE}/esd_list"
}

#==============================
# Copy physics lists
#==============================

copy_run_list() {
    local out_all="${DEST_BASE}/physics_good.txt"
    : > "${out_all}"

    for ph in "${PHASES[@]}"; do
        local f="${ORIG_BASE}/${ph}/goodrunlist_${VERSION}/physics_good.txt"

        if [[ ! -f "${f}" ]]; then
            log WARN "Missing: ${f}"
            continue
        fi

        local out="${DEST_BASE}/physics_good_${ph}.txt"
        cp -f -- "${f}" "${out}"

        cat "${f}" >> "${out_all}"
        log INFO "Merged ${ph}"
    done
}

#==============================
# Copy file list
#==============================

copy_file_list() {

    for ph in "${PHASES[@]}"; do
        local base="${ORIG_BASE}/${ph}/goodrunlist_${VERSION}"

        local rtraw="${base}/rtraw_list"
        local esd="${base}/kup_esd_list"

        if [[ ! -d "${rtraw}" ]]; then
            log WARN "Missing: ${rtraw}"
            continue
        fi
        if [[ ! -d "${esd}" ]]; then
            log WARN "Missing: ${esd}"
            continue
        fi

        local files=( "${rtraw}/run_"*.txt )
        if [[ ! -e "${files[0]}" ]]; then
            log WARN "No run_*.txt files in ${rtraw}"
            continue
        fi
        cp -f -- "${files[@]}" "${DEST_BASE}/rtraw_list/" &

        files=( "${esd}/run_"*.txt )
        if [[ ! -e "${files[0]}" ]]; then
            log WARN "No run_*.txt files in ${esd}"
            continue
        fi
        cp -f -- "${files[@]}" "${DEST_BASE}/esd_list/" &
    done
}

#==============================
# Main
#==============================

main() {
    parse_args "$@"
    
    clean_dir
    copy_run_list
    copy_file_list

    wait
    local status=$?
    if [[ $status -ne 0 ]]; then
        log ERROR "One or more background copies failed (exit: ${status})"
        exit $status
    fi
    
    log INFO "Done"
}

main "$@"