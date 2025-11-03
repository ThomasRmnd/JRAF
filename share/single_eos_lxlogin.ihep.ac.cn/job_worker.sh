#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO Job Worker
#  Purpose: Process a single ESD–RTRAW pair, with optional neighboring files
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

#==============================
# Parse command-line arguments
#==============================

usage() {
    cat <<EOF
Usage: $(basename "$0") <proc_id> <run_number> <list_base> [extra_args...]

Process a single ESD–RTRAW pair identified by run and process ID.

Arguments:
  <proc_id>       Process index within the list
  <run_number>    Run number to process
  <list_base>     Path to the list base directory (contains rtraw_list/esd_list)
  [extra_args...] Optional additional arguments passed to run.py
EOF
}

parse_args() {
    if (( $# < 3 )); then
        log ERROR "Missing required arguments"
        usage >&2
        exit 1
    fi

    PROC_ID="$1"; shift
    RUN_NUMBER="$1"; shift
    LIST_BASE="$1"; shift
    EXTRA_ARGS=("$@")
}

#==============================
# File Preparation
#==============================

load_file_lists() {
    local rtraw_list_file="${LIST_BASE}/rtraw_list/run_${RUN_NUMBER}.txt"
    local esd_list_file="${LIST_BASE}/esd_list/run_${RUN_NUMBER}.txt"

    log INFO "Listing ROOT files from EOS..."
    mapfile -t RTRAW_LIST < <(xrdfs "$EOS_BASE" cat "$rtraw_list_file")
    mapfile -t ESD_LIST   < <(xrdfs "$EOS_BASE" cat "$esd_list_file")

    log INFO "Number of RTRAW files: ${#RTRAW_LIST[@]}"
    log INFO "Number of ESD   files: ${#ESD_LIST[@]}"
}

get_file_number() {
    local fname="$1"
    [[ $fname =~ \.[0-9]{14}\.([0-9]+)_ ]] && echo "${BASH_REMATCH[1]}" || echo ""
}

include_neighbor() {
    local index="$1"
    local direction="$2"
    local step

    case "$direction" in
        prev) step=-1 ;;
        next) step=1 ;;
        *) log WARN "Invalid direction '$direction' in include_neighbor"; return 1 ;;
    esac

    local neighbor=$(( index + step ))

    (( neighbor < 0 || neighbor >= ${#ESD_LIST[@]} )) && return 0

    local fname
    fname=$(basename "${ESD_LIST[$neighbor]}")
    local num
    num=$(get_file_number "$fname")
    num=$((10#$num))

    (( num == target_num + step )) || return 0

    log INFO "Including $direction neighbor: $fname (num=$num)"
    indices_to_process+=("$neighbor")
}

#==============================
# Path Handling
#==============================

resolve_output_paths() {
    local input_esd_file="$1"

    if [[ "$input_esd_file" =~ /eos/juno/esd/([^/]+)/([0-9]{4})/([0-9]{4})/RUN\.([0-9]+)\. ]]; then
        esd_version="${BASH_REMATCH[1]}"
        year="${BASH_REMATCH[2]}"
        monthday="${BASH_REMATCH[3]}"
        output_path="/junofs/users/traymond/analysis/ibd/${year}/${monthday}"
        tt_reco_filepath="${EOS_BASE}/eos/juno/dirac/juno/user/j/jpandre_1/tt_data_auto/${year}/${monthday}/RUN.${RUN_NUMBER}.*.EDM.user.root"
    elif [[ "$input_esd_file" =~ /eos/juno/juno-kup/([^/]+)/([^/]+)/([^/]+)/([^/]+)/([^/]+)/RUN\.([0-9]+)\. ]]; then
        campaign="${BASH_REMATCH[1]}"
        stream="${BASH_REMATCH[2]}"
        run_bucket="${BASH_REMATCH[3]}"
        run_group="${BASH_REMATCH[4]}"
        run_number="${BASH_REMATCH[5]}"
        output_path="/junofs/users/traymond/analysis/ibd/${run_bucket}/${run_group}/${RUN_NUMBER}"
        tt_reco_filepath="${EOS_BASE}/eos/juno/dirac/juno/user/j/jpandre_1/tt_data_auto/${run_bucket}/${run_group}/${RUN_NUMBER}/RUN.${RUN_NUMBER}.*.EDM.user.root"
    else
        log ERROR "Unrecognized ESD path format: $input_esd_file"
        exit 1
    fi

    mkdir -p "$output_path"
    log INFO "Output path: $output_path"
    log DEBUG "TT reco file path: $tt_reco_filepath"
}

#==============================
# Main Execution
#==============================

main() {
    parse_args "$@"
    load_file_lists

    input_esd_file="${ESD_LIST[$PROC_ID]}"
    input_rtraw_file="${RTRAW_LIST[$PROC_ID]}"
    input_esd_filename=$(basename "$input_esd_file")

    target_num=$(get_file_number "$input_esd_filename")
    if [[ -z $target_num ]]; then
        log ERROR "Cannot extract file number from $input_esd_filename"
        exit 1
    fi
    target_num=$((10#$target_num))

    indices_to_process=("$PROC_ID")
    include_neighbor "$PROC_ID" "prev"
    include_neighbor "$PROC_ID" "next"

    indices_to_process=($(printf "%s\n" "${indices_to_process[@]}" | sort -n))
    log INFO "Files to process: ${indices_to_process[*]}"

    resolve_output_paths "$input_esd_file"
    
    source /afs/ihep.ac.cn/users/t/traymond/J25.3.0/git_junosw_J25_load.sh
    log INFO "Environment loaded (TUTORIALROOT=${TUTORIALROOT})"
    log INFO "Temporary directory: ${TEMP}"

    esd_files=()
    rtraw_files=()

    for idx in "${indices_to_process[@]}"; do
        local_esd_file="${TEMP}/$(basename "${ESD_LIST[$idx]}")"
        local_rtraw_file="${TEMP}/$(basename "${RTRAW_LIST[$idx]}")"

        log INFO "Copying input files for index $idx"
        xrdcp "${ESD_LIST[$idx]}" "$local_esd_file"
        xrdcp "${RTRAW_LIST[$idx]}" "$local_rtraw_file"

        esd_files+=("$local_esd_file")
        rtraw_files+=("$local_rtraw_file")
    done

    local_output_file="${TEMP}/${input_esd_filename/.esd/.output.root}"
    output_file="$output_path/$(basename "$local_output_file")"

    log INFO "Running reconstruction with run.py..."
    python run.py \
        --input "${esd_files[@]}" \
        --input-rtraw "${rtraw_files[@]}" \
        --output "$local_output_file" \
        --target-input-filename "$input_esd_filename" \
        --tt-reco-filepath "$tt_reco_filepath" \
        "${EXTRA_ARGS[@]}"

    if cp "$local_output_file" "$output_file"; then
        log INFO "Output copied to $output_file"
    else
        log ERROR "Failed to copy output file to destination"
        exit 1
    fi

    log INFO "Completed run ${RUN_NUMBER} (PROC_ID=${PROC_ID})"
}

main "$@"