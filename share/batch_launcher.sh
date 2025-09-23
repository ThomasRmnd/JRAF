#!/bin/bash
# batch_launcher.sh
#
# Usage:
#   ./batch_launcher.sh -f PhysicsRun.txt
#   ./batch_launcher.sh -f PhysicsRun.txt -r 9737 -R 9754
#
# PhysicsRun.txt layout:
#   run_number output_path

set -euo pipefail

usage() {
    echo "Usage:"
    echo "  $0 -f PhysicsRun.txt"
    echo "  $0 -f PhysicsRun.txt -r <lower> -R <upper>"
    exit 1
}

FILE=""
LOWER=""
UPPER=""

while getopts "f:r:R:" opt; do
    case $opt in
        f) FILE="$OPTARG" ;;
        r) LOWER="$OPTARG" ;;
        R) UPPER="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "$FILE" ]]; then
    usage
fi

while read -r run outpath; do
    [[ -z "$run" || -z "$outpath" ]] && continue

    if [[ -n "$LOWER" && -n "$UPPER" ]]; then
        if (( run < LOWER || run > UPPER )); then
            continue
        fi
    fi

    echo ">>> Launching run $run -> $outpath"
    sh job_launcher.sh --run-number "$run" --output-path "$outpath"

done < "$FILE"
