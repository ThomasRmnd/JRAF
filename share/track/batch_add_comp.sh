#!/bin/bash
#--------------------------------------------------------------------------------------------------
#  JUNO Multi-run Hadd Helper
#  Purpose: hadd all ROOT outputs for multiple runs automatically
#--------------------------------------------------------------------------------------------------

set -euo pipefail
IFS=$'\n\t'

#==============================
# Utility functions
#==============================

source /pbs/home/t/traymond/share/bash/logging.sh

#==============================
# Configuration defaults
#==============================

BASE_OUTPUT_DIR="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/comparison"

LOWER_BOUND=""
UPPER_BOUND=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --lower <num>      Starting run number (inclusive)
  --upper <num>      Ending run number (inclusive)
  --help             Show this help message and exit

Examples:
  $(basename "$0") --lower 9000 --upper 9050
  $(basename "$0") --list /path/to/custom_run_list.txt
EOF
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --lower) LOWER_BOUND="$2"; shift 2 ;;
            --upper) UPPER_BOUND="$2"; shift 2 ;;
            --help|-h) usage; exit 0 ;;
            *) log ERROR "Unknown argument: $1"; usage; exit 1 ;;
        esac
    done

    if [[ -z "$LOWER_BOUND" || -z "$UPPER_BOUND" ]]; then
        log ERROR "Both --lower and --upper must be specified"
        usage
        exit 1
    fi
}

#==============================
# File helpers
#==============================
extract_run_number() {
    local filename="$1"

    if [[ "$filename" =~ ^RUN\.([0-9]+)\.comparison\.output\.reprod25c\.cca\.root$ ]]; then
        echo "${BASH_REMATCH[1]}"
    else
        echo ""
    fi
}

is_run_in_range() {
    local run="$1"
    (( run >= LOWER_BOUND && run <= UPPER_BOUND ))
}

#==============================
# File discovery
#==============================
collect_input_files() {
    local -n _files=$1
    _files=()

    for file in "$BASE_OUTPUT_DIR"/RUN.*.comparison.output.reprod25c.cca.root; do
        [[ -e "$file" ]] || continue

        local filename
        filename=$(basename "$file")

        local run
        run=$(extract_run_number "$filename")
        [[ -z "$run" ]] && continue

        if is_run_in_range "$run"; then
            _files+=("$file")
        fi
    done
}


#==============================
# Hadd execution
#==============================
run_hadd() {
    local -a files=("$@")

    if (( ${#files[@]} == 0 )); then
        log ERROR "No files found in requested run range"
        exit 1
    fi

    local output_file
    output_file="$BASE_OUTPUT_DIR/RUN.${LOWER_BOUND}-${UPPER_BOUND}.comparison.output.reprod25c.cca.root"

    log INFO "Found ${#files[@]} files"
    log INFO "Output file: $output_file"
    log INFO "Running hadd..."

    hadd -f "$output_file" "${files[@]}"

    log INFO "Hadd completed successfully"
}

#==============================
# Main
#==============================
main() {
    parse_args "$@"
    validate_inputs

    log INFO "Searching runs in range [$LOWER_BOUND, $UPPER_BOUND]"
    log INFO "Directory: $BASE_OUTPUT_DIR"

    local -a input_files
    collect_input_files input_files

    run_hadd "${input_files[@]}"
}

main "$@"