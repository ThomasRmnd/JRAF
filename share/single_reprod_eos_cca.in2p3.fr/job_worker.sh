#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO Job Worker
#  Purpose: Process a single ReProd, with optional neighboring files
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

EOS_BASE="root://junoeos01.ihep.ac.cn/"

#==============================
# Parse command-line arguments
#==============================

usage() {
    cat <<EOF
Usage: $(basename "$0") <run_number> <list_base> [extra_args...]

Process a single ReProd identified by run and process ID.

Arguments:
  <run_number>    Run number to process
  <list_base>     Path to the list base directory (contains rtraw_list)
  [extra_args...] Optional additional arguments passed to run.py
EOF
}

parse_args() {
    if (( $# < 3 )); then
        log ERROR "Missing required arguments"
        usage >&2
        exit 1
    fi

    PROC_ID=${SLURM_ARRAY_TASK_ID}
    RUN_NUMBER="$1"; shift
    LIST_BASE="$1"; shift
    EXTRA_ARGS=("$@")
}

#==============================
# File Preparation
#==============================

load_file_lists() {
    local rtraw_list_file="${LIST_BASE}/rtraw_list/run_${RUN_NUMBER}.txt"

    log INFO "Listing ROOT files from EOS..."
    mapfile -t RTRAW_LIST   < <(xrdfs "$EOS_BASE" cat "$rtraw_list_file")

    log INFO "Number of RTRAW files: ${#RTRAW_LIST[@]}"
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

    (( neighbor < 0 || neighbor >= ${#RTRAW_LIST[@]} )) && return 0

    local fname
    fname=$(basename "${RTRAW_LIST[$neighbor]}")
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
    local input_reprod_file="$1"

    if [[ "$input_reprod_file" =~ /eos/juno/juno-reprod/([^/]+)/([^/]+)/([^/]+)/([^/]+)/([^/]+)/RUN\.([0-9]+)\. ]]; then
        campaign="${BASH_REMATCH[1]}"
        stream="${BASH_REMATCH[2]}"
        run_bucket="${BASH_REMATCH[3]}"
        run_group="${BASH_REMATCH[4]}"
        run_number="${BASH_REMATCH[5]}"
        output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/ibd/${run_bucket}/${run_group}/${RUN_NUMBER}"
        # tt_reco_filepath="${EOS_BASE}/eos/juno/dirac/juno/user/j/jpandre_1/tt_data_auto/${run_bucket}/${run_group}/${RUN_NUMBER}/RUN.${RUN_NUMBER}.*.EDM.user.root"
        tt_reco_filepath=""
    else
        log ERROR "Unrecognized ReProd path format: $input_reprod_file"
        exit 1
    fi

    mkdir -p "$output_path"
    log INFO "Output path: $output_path"
}

#==============================
# RTRAW to ReProd filename
#==============================

rtraw_to_reprod_filename() {
    local fpath="$1"
    local fname
    fname=$(basename "$fpath")

    local run stream run_bucket run_group

    if [[ "$fpath" =~ /eos/juno/rtraw/([0-9]{4})/([0-9]{4})/RUN\.([0-9]+)\. ]]; then
        run="${BASH_REMATCH[3]}"
        bucket_val=$(( (10#$run / 1000) * 1000 ))
        group_val=$(( (10#$run / 100) * 100 ))
        run_bucket=$(printf "%08d" "$bucket_val")
        run_group=$(printf "%08d" "$group_val")

    elif [[ "$fpath" =~ /eos/juno/juno-rtraw/([^/]+)/([^/]+)/([0-9]+)/([0-9]+)/([0-9]+)/RUN\.([0-9]+)\. ]]; then
        run="${BASH_REMATCH[6]}"
        stream="${BASH_REMATCH[2]}"
        run_bucket="${BASH_REMATCH[3]}"
        run_group="${BASH_REMATCH[4]}"
    else
        echo "ERROR: can't extract run number or structure from path: $fpath" >&2
        exit 1
    fi

    if [[ -z "${stream:-}" ]]; then
        if [[ "$fname" =~ \.([^.]+)\.[0-9]{14}\.[0-9]+_ ]]; then
            stream="${BASH_REMATCH[1]}"
        else
            stream="unknown_stream"
        fi
    fi

    local output_reprod_filename="${fname/.rtraw/.esd}"

    local base_dir="/eos/juno/juno-reprod/ReProd25B/${stream}/${run_bucket}"
    local eos_base="root://junoeos01.ihep.ac.cn/"

    local candidate_groups
    candidate_groups=$(xrdfs "$eos_base" ls "$base_dir" 2>/dev/null | grep -E "/${run_group}(_v[0-9]+)?/?$" | sort)

    if [[ -z "$candidate_groups" ]]; then
        echo "ERROR: no run_group directory found under $base_dir" >&2
        exit 1
    fi

    local selected_group
    selected_group=$(basename "$(echo "$candidate_groups" | tail -n 1)")

    local reprod_path="${eos_base}/eos/juno/juno-reprod/ReProd25B/${stream}/${run_bucket}/${selected_group}/${run}/${output_reprod_filename}"

    echo "$reprod_path"
}


#==============================
# Main Execution
#==============================

main() {
    parse_args "$@"
    load_file_lists

    input_rtraw_file="${RTRAW_LIST[$PROC_ID]}"
    input_rtraw_filename=$(basename "$input_rtraw_file")

    target_num=$(get_file_number "$input_rtraw_filename")
    if [[ -z $target_num ]]; then
        log ERROR "Cannot extract file number from $input_rtraw_filename"
        exit 1
    fi
    target_num=$((10#$target_num))

    indices_to_process=("$PROC_ID")
    include_neighbor "$PROC_ID" "prev"
    include_neighbor "$PROC_ID" "next"

    indices_to_process=($(printf "%s\n" "${indices_to_process[@]}" | sort -n))
    log INFO "Files to process: ${indices_to_process[*]}"

    input_reprod_file=$(rtraw_to_reprod_filename "$input_rtraw_file")
    
    source /pbs/home/t/traymond/J25.6.1_Modified/git_junosw_J25_load.sh
    log INFO "Environment loaded (TUTORIALROOT=${TUTORIALROOT})"
    log INFO "Temporary directory: ${TMPDIR}"

    reprod_files=()

    for idx in "${indices_to_process[@]}"; do
        input_rtraw_file="${RTRAW_LIST[$idx]}"
        input_reprod_file=$(rtraw_to_reprod_filename "$input_rtraw_file")
        local_reprod_file="${TMPDIR}/$(basename "${input_reprod_file}")"

        log INFO "Copying input files for index $idx"
        xrdcp "${input_reprod_file}" "$local_reprod_file"

        reprod_files+=("$local_reprod_file")
    done

    input_rtraw_file="${RTRAW_LIST[$PROC_ID]}"
    input_rtraw_filename=$(basename "$input_rtraw_file")
    input_reprod_file=$(rtraw_to_reprod_filename "$input_rtraw_file")
    input_reprod_filename=$(basename "$input_reprod_file")

    log INFO "Input file: $input_reprod_file"

    resolve_output_paths "$input_reprod_file"

    local_output_file="${TMPDIR}/${input_reprod_filename/.esd/.output.root}"
    output_file="$output_path/$(basename "$local_output_file")"

    log INFO "Running reconstruction with run.py..."
    python run.py \
        --input "${reprod_files[@]}" \
        --output "$local_output_file" \
        --target-input-filename "$input_reprod_filename" \
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