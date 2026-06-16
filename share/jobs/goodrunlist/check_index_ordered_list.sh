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
# Configuration defaults
#==============================

GOODRUNLIST_PATH="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/other/GoodList"

#==============================
# Usage & Argument Parsing
#==============================

VERSION=""
CAMPAIGN=""

usage() {
    cat <<EOF
Usage: $(basename "$0") --campaign <str> --version <str>

Required:
  --campaign            <str>           Campaign selection {ReProd25C|ReProd25D|ValProd26B|ReProd26B}
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
            --campaign)     CAMPAIGN="$2"; shift 2 ;;
            --version)       VERSION="$2"; shift 2 ;;
            --help|-h)       usage; exit 0 ;;
            *) log ERROR "Unknown argument: $1"; usage; exit 1 ;;
        esac
    done

    if [[ -z "${CAMPAIGN:-}" ]]; then
        log ERROR "--campaign is required {ReProd25C|ReProd25D|ValProd26B|ReProd26B}"
        usage
        exit 1
    fi

    case "${CAMPAIGN}" in
        ReProd25C|ReProd25D|ValProd26B|ReProd26B) ;;
        *) log ERROR "Invalid --campaign: ${CAMPAIGN} (expected {ReProd25C|ReProd25D|ValProd26B|ReProd26B})"
           exit 1 ;;
    esac

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
# Check
#==============================

check_file_order() {
    local filepath="$1"
    local filename
    filename="$(basename "${filepath}")"
 
    local index_re='\.([0-9]+)_J[0-9]'
 
    local prev_index=-1
    local line_num=0
    local ok=0
 
    while IFS= read -r line || [[ -n "${line}" ]]; do
        (( line_num++ )) || true
 
        if [[ -z "${line}" ]]; then
            continue
        fi
 
        if [[ "${line}" =~ ${index_re} ]]; then
            local raw_index="${BASH_REMATCH[1]}"
            local index=$(( 10#${raw_index} ))
 
            if (( index <= prev_index )); then
                log WARN "${filename}: ordering violation at line ${line_num}: index ${raw_index} follows ${prev_index}"
                ok=1
            fi
 
            prev_index="${index}"
        else
            log WARN "${filename}: line ${line_num} has no recognisable file index — skipped"
            log DEBUG "  -> ${line}"
        fi
    done < "${filepath}"
 
    if (( line_num == 0 )); then
        log WARN "${filename}: file is empty"
        return 1
    fi
 
    return "${ok}"
}

check_list_dir() {
    local dir="$1"
    local label="$2"
    local issues=0
 
    if [[ ! -d "${dir}" ]]; then
        log WARN "Directory not found: ${dir}"
        return 1
    fi
 
    local run_files=( "${dir}/run_"*.txt )
 
    if [[ ! -e "${run_files[0]}" ]]; then
        log WARN "No run_*.txt files in: ${dir}"
        return 1
    fi
 
    local total="${#run_files[@]}"
    local checked=0
 
    for f in "${run_files[@]}"; do
        if [[ ! -f "${f}" ]]; then
            log WARN "Not a regular file, skipping: ${f}"
            continue
        fi
 
        (( checked++ )) || true
 
        if ! check_file_order "${f}"; then
            (( issues++ )) || true
        fi
    done
 
    if (( issues == 0 )); then
        log INFO "${label}: all ${checked}/${total} files are correctly ordered"
    else
        log WARN "${label}: ${issues}/${checked} file(s) have ordering issues"
    fi
 
    return $(( issues > 0 ? 1 : 0 ))
}

run_checks() {
    local global_rc=0
    local base="${GOODRUNLIST_PATH}/${CAMPAIGN}"

    if [[ ! -d "${base}" ]]; then
        log WARN "Phase directory not found: ${base}"
        return 1
    fi

    local rtraw_dir="${base}/rtraw_list"
    check_list_dir "${rtraw_dir}" "rtraw_list" || global_rc=1

    local esd_dir="${base}/esd_list"
    check_list_dir "${esd_dir}" "esd_list" || global_rc=1
 
    return "${global_rc}"
}

#==============================
# Main
#==============================

main() {
    parse_args "$@"
 
    local rc=0
    run_checks || rc=$?
 
    if (( rc == 0 )); then
        log INFO "All file lists are correctly ordered"
    else
        log WARN "One or more file lists have ordering issues. See warnings above"
    fi
 
    exit "${rc}"
}

main "$@"