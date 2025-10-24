#!/bin/bash
# batch_launcher.sh
#
# Usage:
#   ./batch_launcher.sh
#   ./batch_launcher.sh -r <num> -R <num>

set -euo pipefail

EOS_BASE="root://junoeos01.ihep.ac.cn/"
LIST_BASE="/eos/juno/groups/DataQuality/P25A/Physics/goodrunlist_v3.4/Physics_good_run_list.txt"

usage() {
    echo "Usage:"
    echo "  $0 [-r <lower>] [-R <upper>]"
    exit 1
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

echo "Fetching run list from EOS..."
mapfile -t RUN_LIST < <(xrdfs "$EOS_BASE" cat "$LIST_BASE" | tr -d '\r' | sed '/^$/d')

echo "Total runs found: ${#RUN_LIST[@]}"

for run in "${RUN_LIST[@]}"; do
    # Ensure run is numeric (ignore headers or malformed lines)
    if ! [[ "$run" =~ ^[0-9]+$ ]]; then
        continue
    fi

    # Apply optional range filter
    if [[ -n "$LOWER" && -n "$UPPER" ]]; then
        if (( run < LOWER || run > UPPER )); then
            continue
        fi
    fi

    echo ">>> Launching run $run"
    sh job_launcher.sh --run-number "$run"
done
