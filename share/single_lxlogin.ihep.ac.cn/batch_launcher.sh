#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO job multi-submission helper
#  Purpose: Automate hep_sub job multi-submissions for single ESD–RTRAW processing
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
LIST_BASE="/eos/juno/groups/DataQuality/P25A/Physics/goodrunlist_v3.4/Physics_good_run_list.txt"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Optional:
  --r <num>                    Starting run number (default: the first in the run list)
  --R <num>                    Ending run number (default: the last in the run list)
  --help                       Show this help message and exit
EOF
}

LOWER=""
UPPER=""

while getopts "r:R:" opt; do
    case $opt in
        r) LOWER="$OPTARG" ;;
        R) UPPER="$OPTARG" ;;
        *) usage ;;
    esac
done

log INFO "Fetching run list from EOS..."
mapfile -t RUN_LIST < <(xrdfs "$EOS_BASE" cat "$LIST_BASE" | tr -d '\r' | sed '/^$/d')

log INFO "Total runs found: ${#RUN_LIST[@]}"

for run in "${RUN_LIST[@]}"; do
    if ! [[ "$run" =~ ^[0-9]+$ ]]; then
        continue
    fi

    if [[ -n "$LOWER" && -n "$UPPER" ]]; then
        if (( run < LOWER || run > UPPER )); then
            continue
        fi
    fi

    log INFO ">>> Launching run $run"
    sh job_launcher.sh --run-number "$run"
done
