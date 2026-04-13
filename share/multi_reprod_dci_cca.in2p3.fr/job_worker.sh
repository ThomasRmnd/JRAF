#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO Job Worker
#  Purpose: Process multi ReProd, with optional neighboring files
#--------------------------------------------------------------------------------------------------

set -euo pipefail
IFS=$'\n\t'

#==============================
# Utility functions
#==============================

source /pbs/home/t/traymond/share/bash/logging.sh
export X509_USER_PROXY=/sps/juno/jdeandre/rtraw_ThomasRaymond/.cert_traymond_juno_user

#==============================
# Configuration defaults
#==============================

XRD_URL_EOS="root://junoeos01.ihep.ac.cn/"
XRD_URL_CNAF="root://xrootd-archive.cr.cnaf.infn.it:1095/"

XRD_BASEPATH_EOS="/eos"
XRD_BASEPATH_CNAF="/production/storm/dirac"
LOCAL_BASEPATH="/sps/juno/jdeandre/rtraw_ThomasRaymond"

OUTPUT_SUFFIX_NORMAL="output.normal.root"
OUTPUT_SUFFIX_REPROD25A="output.reprod25a.root"
OUTPUT_SUFFIX_REPROD25B="output.reprod25b.root"
OUTPUT_SUFFIX_REPROD25C="output.reprod25c.root"
OUTPUT_SUFFIX_REPROD25D="output.reprod25d.root"

#==============================
# Global Flags
#==============================

DIRECT_IO=0 # default = copy to TMPDIR
SKIP_IF_EXIST=0 # default = doesn't skip
USE_LOCAL=0 # default = use remote xrootd

#==============================
# Parse command-line arguments
#==============================

usage() {
    cat <<EOF
Usage: $(basename "$0") <site> <campaign> <run> <list-base> <range-start> <range-end> [options...]

Process a single range of files identified by run and file indices.

Arguments:
  <str>                          Storage site selection {EOS|CNAF}
  <str>                          Campaign selection {Normal|ReProd25A|ReProd25B|ReProd25C|ReProd25D}
  <int>                          Run number to process
  <str>                          Basepath for the file list
  <int>                          Start index of the file range in the run list
  <int>                          End index of the file range in the run list

Options:
  --local                        Use local files instead of remote xrootd
  --no-local-copy | --direct-io  Use direct I/O (no copy to TMPDIR)
  --skip-if-exist                Skip the job if the final output file already exists.
EOF
}

parse_args() {
    if (( $# < 5 )); then
        log ERROR "Missing required arguments"
        usage >&2
        exit 1
    fi

    SITE="$1"; shift
    CAMPAIGN="$1"; shift
    RUN_NUMBER="$1"; shift
    LIST_BASE="$1"; shift
    RANGE_START="$1"; shift
    RANGE_END="$1"; shift

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --local)
                USE_LOCAL=1
                ;;
            --no-local-copy|--direct-io)
                DIRECT_IO=1
                ;;
            --skip-if-exist)
                SKIP_IF_EXIST=1
                ;;
            *)
                EXTRA_ARGS+=("$1")
                ;;
        esac
        shift
    done

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
            OUTPUT_SUFFIX="${OUTPUT_SUFFIX_NORMAL}"
            ;;
        ReProd25A)
            OUTPUT_SUFFIX="${OUTPUT_SUFFIX_REPROD25A}"
            ;;
        ReProd25B)
            OUTPUT_SUFFIX="${OUTPUT_SUFFIX_REPROD25B}"
            ;;
        ReProd25C)
            OUTPUT_SUFFIX="${OUTPUT_SUFFIX_REPROD25C}"
            ;;
        ReProd25D)
            OUTPUT_SUFFIX="${OUTPUT_SUFFIX_REPROD25D}"
            ;;
        *)
            log ERROR "Invalid campaign argument: ${CAMPAIGN} (expected {Normal|ReProd25A|ReProd25B|ReProd25C})"
            exit 1
            ;;
    esac
}

#==============================
# File Preparation
#==============================

load_file_lists() {
    local rtraw_list_file="${LIST_BASE}/rtraw_list/run_${RUN_NUMBER}.txt"

    log INFO "Listing ROOT files from EOS..."
    mapfile -t RTRAW_LIST   < <(cat "${rtraw_list_file}")

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

    local fname_this fname_neighbor num_this num_neighbor

    fname_this=$(basename "${RTRAW_LIST[$index]}")
    fname_neighbor=$(basename "${RTRAW_LIST[$neighbor]}")

    num_this=$(get_file_number "$fname_this")
    num_neighbor=$(get_file_number "$fname_neighbor")

    if [[ -z "$num_this" || -z "$num_neighbor" ]]; then
        log ERROR "Neighbor check: couldn't parse numbers for index $index or $neighbor: '$fname_this' / '$fname_neighbor'"
        return 1
    fi

    num_this=$((10#$num_this))
    num_neighbor=$((10#$num_neighbor))

    (( num_neighbor == num_this + step )) || return 0

    log INFO "Including $direction neighbor: $fname_neighbor (num=$num_neighbor)"
    indices_to_process+=("$neighbor")
}


#==============================
# Path Handling
#==============================

resolve_input_paths() {
    local input_reprod_file="$1"

    if [[ "$input_reprod_file" =~ /juno/juno-reprod/([^/]+)/([^/]+)/([^/]+)/([^/]+)/([^/]+)/RUN\.([0-9]+)\.[^/]*([0-9]{14})[^/]* ]]; then
        campaign="${BASH_REMATCH[1]}"
        stream="${BASH_REMATCH[2]}"
        run_bucket="${BASH_REMATCH[3]}"
        run_group="${BASH_REMATCH[4]}"
        run_number="${BASH_REMATCH[5]}"
        timestamp="${BASH_REMATCH[7]}"
    else
        log ERROR "Unrecognized ReProd path format: $input_reprod_file"
        exit 1
    fi

    if (( RUN_NUMBER >= 9591 && RUN_NUMBER <= 10169 )); then
        year="${timestamp:0:4}"
        month="${timestamp:4:2}"
        day="${timestamp:6:2}"
        tt_reco_filepath="${XRD_URL}${XRD_BASEPATH}/juno/user/j/jpandre_1/tt_data_auto/${year}/${month}${day}/RUN.${RUN_NUMBER}.*.EDM.user.root"

    elif (( RUN_NUMBER >= 10176 && RUN_NUMBER <= 10479 )); then
        tt_reco_filepath="${XRD_URL}${XRD_BASEPATH}/juno/user/j/jpandre_1/tt_data_auto/${run_bucket}/${run_group}/${RUN_NUMBER}/RUN.${RUN_NUMBER}.*.EDM.user.root"

    elif (( RUN_NUMBER >= 10480 )); then
        tt_reco_filepath="${XRD_URL}${XRD_BASEPATH}/juno/juno-reprod/TT25A/J25.4.3-patched/user_rec/${run_bucket}/${run_group}/${RUN_NUMBER}/RUN.${RUN_NUMBER}.*.EDM.user.root"

    else
        log ERROR "No TT reco path rule defined for run ${RUN_NUMBER}"
        exit 1
    fi

    log INFO "TT-Reco filepath resolved: ${tt_reco_filepath}"
}

resolve_output_paths() {
    local input_reprod_file="$1"

    if [[ "$input_reprod_file" =~ /juno/juno-reprod/([^/]+)/([^/]+)/([^/]+)/([^/]+)/([^/]+)/RUN\.([0-9]+)\. ]]; then
        campaign="${BASH_REMATCH[1]}"
        stream="${BASH_REMATCH[2]}"
        run_bucket="${BASH_REMATCH[3]}"
        run_group="${BASH_REMATCH[4]}"
        run_number="${BASH_REMATCH[5]}"
        output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/ibd/${run_bucket}/${run_group}/${RUN_NUMBER}"
        reco_output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/${run_bucket}/${run_group}/${RUN_NUMBER}"
        # feature_output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/features/reprod/${run_bucket}/${run_group}/${RUN_NUMBER}"
    else
        log ERROR "Unrecognized ReProd path format: $input_reprod_file"
        exit 1
    fi

    mkdir -p "${output_path}"
    mkdir -p "${reco_output_path}"
    # mkdir -p "${feature_output_path}"

    log INFO "Output path: ${output_path}"
    log DEBUG "TT reco file path: ${tt_reco_filepath}"
}

check_output_existence() {
    local input_rtraw_file="${RTRAW_LIST[$RANGE_START]}"
    local input_reprod_file=$(rtraw_to_reprod_filename "${input_rtraw_file}")

    local input_reprod_filename=$(basename "${input_reprod_file}")

    if ! resolve_output_paths "${input_reprod_file}"; then
        log ERROR "Failed to resolve output paths for existence check"
        exit 1
    fi
    
    local output_filename="RUN.${RUN_NUMBER}.${RANGE_START}-${RANGE_END}.${OUTPUT_SUFFIX}"
    local output_file="${output_path}/${output_filename}"
    
    if [[ -f "${output_file}" ]]; then
        log INFO "Output file already exists, skipping job: ${output_file}"
        exit 0
    fi
}

#==============================
# RTRAW to ReProd filename
#==============================

rtraw_to_reprod_filename() {
    local fpath="$1"
    local fname
    fname=$(basename "$fpath")

    local run stream run_bucket run_group

    if [[ "$fpath" =~ /juno/rtraw/([0-9]{4})/([0-9]{4})/RUN\.([0-9]+)\. ]]; then
        run="${BASH_REMATCH[3]}"
        bucket_val=$(( (10#$run / 1000) * 1000 ))
        group_val=$(( (10#$run / 100) * 100 ))
        run_bucket=$(printf "%08d" "$bucket_val")
        run_group=$(printf "%08d" "$group_val")

    elif [[ "$fpath" =~ /juno/juno-rtraw/([^/]+)/([^/]+)/([0-9]+)/([0-9]+)/([0-9]+)/RUN\.([0-9]+)\. ]]; then
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

    local base_dir="${XRD_BASEPATH}/juno/juno-reprod/${CAMPAIGN}/${stream}/${run_bucket}"

    local candidate_groups
    candidate_groups=$(xrdfs "${XRD_URL}" ls "${base_dir}" 2>/dev/null | grep -E "/${run_group}(_v[0-9]+)?/?$" | sort)

    if [[ -z "${candidate_groups}" ]]; then
        log ERROR "No run_group (${run_group}) directory found under ${base_dir}" >&2
        exit 1
    fi

    local selected_group
    selected_group=$(basename "$(echo "$candidate_groups" | tail -n 1)")

    local reprod_path=""
    if (( USE_LOCAL == 1 )); then
        reprod_path="${LOCAL_BASEPATH}/juno/juno-reprod/${CAMPAIGN}/${stream}/${run_bucket}/${selected_group}/${run}/${output_reprod_filename}"
    else
        reprod_path="${XRD_URL}${XRD_BASEPATH}/juno/juno-reprod/${CAMPAIGN}/${stream}/${run_bucket}/${selected_group}/${run}/${output_reprod_filename}"
    fi

    echo "${reprod_path}"
}


#==============================
# Main Execution
#==============================

main() {
    parse_args "$@"
    load_file_lists

    if (( SKIP_IF_EXIST == 1 )); then
        check_output_existence
    fi

    input_first_rtraw_file="${RTRAW_LIST[$RANGE_START]}"
    input_first_rtraw_filename=$(basename "${input_first_rtraw_file}")
    input_last_rtraw_file="${RTRAW_LIST[$RANGE_END]}"
    input_last_rtraw_filename=$(basename "${input_last_rtraw_file}")

    indices_to_process=()

    include_neighbor "${RANGE_START}" "prev"
    prev_idx=$((RANGE_START - 1))
    if (( prev_idx >= 0 )); then
        prev_file_local="$(basename "$(rtraw_to_reprod_filename "${RTRAW_LIST[$prev_idx]}")")"
    else
        prev_file_local=""
    fi

    for ((k = RANGE_START; k <= RANGE_END; k++)); do
        indices_to_process+=("$k")
    done

    include_neighbor "${RANGE_END}" "next"
    next_idx=$((RANGE_END + 1))
    if (( next_idx < ${#RTRAW_LIST[@]} )); then
        next_file_local="$(basename "$(rtraw_to_reprod_filename "${RTRAW_LIST[$next_idx]}")")"
    else
        next_file_local=""
    fi

    indices_to_process=($(printf "%s\n" "${indices_to_process[@]}" | sort -n))
    log INFO "Files to process: ${indices_to_process[*]}"

    source /pbs/home/t/traymond/J25.7.4/git_junosw_load_J25_7_4.sh
    log INFO "Environment loaded (TUTORIALROOT=${TUTORIALROOT})"
    log INFO "Temporary directory: ${TMPDIR}"

    reprod_files=()

    for idx in "${indices_to_process[@]}"; do
        input_rtraw_file="${RTRAW_LIST[$idx]}"
        input_reprod_file=$(rtraw_to_reprod_filename "${input_rtraw_file}")
        filename="$(basename "${input_reprod_file}")"
        local_reprod_file="${TMPDIR}/${filename}"

        if (( DIRECT_IO == 1 )); then
            log INFO "[DIRECT-IO] Using remote file: $filename"
            reprod_files+=("${input_reprod_file}")

        else
            log INFO "Copying input file index ${idx}: $filename"
            if (( USE_LOCAL == 1 )); then
                cp "${input_reprod_file}" "${local_reprod_file}"
            else
                xrdcp "${input_reprod_file}" "${local_reprod_file}"
            fi
            reprod_files+=("${local_reprod_file}")
        fi
    done

    input_rtraw_file="${RTRAW_LIST[$RANGE_START]}"
    input_reprod_file=$(rtraw_to_reprod_filename "${input_rtraw_file}")

    resolve_input_paths "${input_reprod_file}"
    resolve_output_paths "${input_reprod_file}"

    local_output_file="${TMPDIR}/RUN.${RUN_NUMBER}.${RANGE_START}-${RANGE_END}.output.root"
    output_file="$output_path/$(basename "${local_output_file}")"

    local_reco_output_file="${TMPDIR}/RUN.${RUN_NUMBER}.${RANGE_START}-${RANGE_END}.reco.output.root"
    reco_output_file="${reco_output_path}/$(basename "${local_reco_output_file}")"

    # local_feature_output_file="${TMPDIR}/RUN.${RUN_NUMBER}.${RANGE_START}-${RANGE_END}.feature.output.root"
    # feature_output_file="${feature_output_path}/$(basename "${local_feature_output_file}")"

    log INFO "Context previous file: ${prev_file_local:-<none>}"
    log INFO "Context next file: ${next_file_local:-<none>}"

    log INFO "Running run.py..."
    python run.py \
        --input "${reprod_files[@]}" \
        --output "${local_output_file}" \
        --context-previous-filename "${prev_file_local}" \
        --context-next-filename "${next_file_local}" \
        --tt-reco-filepath "${tt_reco_filepath}" \
        --reco-output "${local_reco_output_file}" \
        "${EXTRA_ARGS[@]}"
        # --tt-reco-filepath "${tt_reco_filepath}" \
        # --feature-output "${local_feature_output_file}" \

    if cp "${local_output_file}" "${output_file}"; then
        log INFO "Output copied to ${output_file}"
    else
        log ERROR "Failed to copy ${local_output_file} to ${output_file}"
        exit 1
    fi

    if cp "${local_reco_output_file}" "${reco_output_file}"; then
        log INFO "Output copied to ${reco_output_file}"
    else
        log ERROR "Failed to copy ${local_reco_output_file} to ${reco_output_file}"
        exit 1
    fi

    # if cp "${local_feature_output_file}" "${feature_output_file}"; then
    #     log INFO "Output copied to ${feature_output_file}"
    # else
    #     log ERROR "Failed to copy ${local_feature_output_file} to ${feature_output_file}"
    #     exit 1
    # fi

    log INFO "Completed run ${RUN_NUMBER} (RANGE=[${RANGE_START}, ${RANGE_END}])"
}

main "$@"