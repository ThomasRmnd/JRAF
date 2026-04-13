#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO Job Worker
#  Purpose: Process a single file, with optional neighboring files
#--------------------------------------------------------------------------------------------------

set -euo pipefail
IFS=$'\n\t'

#==============================
# Utility functions
#==============================

HOSTNAME=$(hostname -f 2>/dev/null || hostname)
if [[ "${HOSTNAME}" =~ ^cc.*\.in2p3\.fr$ ]]; then
    # Detect CC-IN2P3 cluster
    CLUSTER="CC-IN2P3"
    source /pbs/home/t/traymond/share/bash/logging.sh
    TEMPDIR=${TMPDIR}
    SOURCE_JUNOSW_PATH="/pbs/home/t/traymond/J25.7.4/git_junosw_load_J25_7_4.sh"
elif [[ "${HOSTNAME}" =~ ^lxlogin[0-9]+\.ihep\.ac\.cn$ ]]; then
    # Detect IHEP cluster
    CLUSTER="IHEP"
    source /junofs/users/traymond/bash/logging.sh
    TEMPDIR=${TEMP}
    SOURCE_JUNOSW_PATH="/afs/ihep.ac.cn/users/t/traymond/J25.3.0/git_junosw_J25_load.sh"
else
    echo "ERROR: Unknown cluster. Hostname: ${HOSTNAME}" >&2
    echo "Expected CC-IN2P3 (cca###) or IHEP (lxlogin###.ihep.ac.cn)" >&2
    exit 1
fi

log INFO "Cluster detected: ${CLUSTER}"

#==============================
# Configuration defaults
#==============================

XRD_URL_EOS="root://junoeos01.ihep.ac.cn/"
XRD_URL_CNAF="root://xrootd-archive.cr.cnaf.infn.it:1095/"

XRD_BASEPATH_EOS="/eos"
XRD_BASEPATH_CNAF="/production/storm/dirac"

LIST_BASE="/eos/juno/groups/DataQuality/P25A/Physics/goodrunlist_v3.6"

#==============================
# Parse command-line arguments
#==============================

usage() {
    cat <<EOF
Usage:
  CC-IN2P3 (SLURM):
      $(basename "$0") <run_number> <site> <campaign> [extra_args...]

  IHEP (hep_sub):
      $(basename "$0") <proc_id> <run_number> <site> <campaign> [extra_args...]

Process a single file identified by run and process ID.

Arguments:
  <proc_id>       Process index (IHEP only)
  <run_number>    Run number to process
  <site>          Storage site selection {EOS|CNAF}
  <campaign>      Campaign selection {Normal|ReProd25A|ReProd25B|ReProd25C}
  [extra_args...] Optional additional arguments passed to run.py
EOF
}

parse_args() {
    if [[ "${CLUSTER}" == "CC-IN2P3" ]]; then
        if (( $# < 3 )); then
            echo "ERROR: Missing arguments for CC-IN2P3" >&2
            usage >&2
            exit 1
        fi
        PROC_ID="${SLURM_ARRAY_TASK_ID}"
    elif [[ "${CLUSTER}" == "IHEP" ]]; then
        if (( $# < 4 )); then
            echo "ERROR: Missing arguments for IHEP" >&2
            usage >&2
            exit 1
        fi
        PROC_ID="$1"; shift
    else
        echo "ERROR: Unknown cluster '${CLUSTER}'" >&2
        exit 1
    fi

    RUN_NUMBER="$1"; shift
    SITE="$1"; shift
    CAMPAIGN="$1"; shift
    EXTRA_ARGS=("$@")

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
            USE_CORRELATION=true
            ;;
        ReProd25A|ReProd25B|ReProd25C)
            USE_CORRELATION=false
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
    local esd_list_file="${LIST_BASE}/esd_list/run_${RUN_NUMBER}.txt"

    log INFO "Listing ROOT files from EOS..."
    mapfile -t RTRAW_LIST < <(xrdfs "${XRD_URL_EOS}" cat "${rtraw_list_file}")
    mapfile -t ESD_LIST   < <(xrdfs "${XRD_URL_EOS}" cat "${esd_list_file}")

    log INFO "Number of RTRAW files: ${#RTRAW_LIST[@]}"
    log INFO "Number of ESD   files: ${#ESD_LIST[@]}"
}

get_file_number() {
    local fname="$1"
    [[ ${fname} =~ \.[0-9]{14}\.([0-9]+)_ ]] && echo "${BASH_REMATCH[1]}" || echo ""
}

include_neighbor() {
    local index="$1"
    local direction="$2"
    local step

    case "${direction}" in
        prev) step=-1 ;;
        next) step=1 ;;
        *) log WARN "Invalid direction '${direction}' in include_neighbor"; return 1 ;;
    esac

    local neighbor=$(( index + step ))

    (( neighbor < 0 || neighbor >= ${#RTRAW_LIST[@]} )) && return 0

    local fname
    fname=$(basename "${RTRAW_LIST[${neighbor}]}")
    local num
    num=$(get_file_number "${fname}")
    num=$((10#$num))

    (( num == target_num + step )) || return 0

    log INFO "Including ${direction} neighbor: ${fname} (num=${num})"
    indices_to_process+=("${neighbor}")
}

#==============================
# Path Handling
#==============================

resolve_output_paths() {
    local input_file="$1"

    if [[ "${input_file}" =~ /juno/rtraw/([0-9]{4})/([0-9]{4})/RUN\.${RUN_NUMBER}\. ]]; then
        local year="${BASH_REMATCH[1]}"
        local monthday="${BASH_REMATCH[2]}"
        local bucket_val=$(( (10#$RUN_NUMBER / 1000) * 1000 ))
        local group_val=$(( (10#$RUN_NUMBER / 100) * 100 ))
        local run_bucket=$(printf "%08d" "$bucket_val")
        local run_group=$(printf "%08d" "$group_val")
        output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/ibd/${run_bucket}/${run_group}/${RUN_NUMBER}"
        reco_output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/${run_bucket}/${run_group}/${RUN_NUMBER}"
        feature_output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/features/reprod/${run_bucket}/${run_group}/${RUN_NUMBER}"
        tt_reco_filepath="${XRD_URL}${XRD_BASEPATH}/juno/user/j/jpandre_1/tt_data_auto/${year}/${monthday}/RUN.${RUN_NUMBER}.*.EDM.user.root"
    elif [[ "${input_file}" =~ /juno/juno-([^/]+)/([^/]+)/([^/]+)/([^/]+)/([^/]+)/${RUN_NUMBER}/RUN\.${RUN_NUMBER}\. ]]; then
        local type="${BASH_REMATCH[1]}"
        local campaign="${BASH_REMATCH[2]}"
        local stream="${BASH_REMATCH[3]}"
        local run_bucket="${BASH_REMATCH[4]}"
        local run_group="${BASH_REMATCH[5]}"
        output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/ibd/${run_bucket}/${run_group}/${RUN_NUMBER}"
        reco_output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/${run_bucket}/${run_group}/${RUN_NUMBER}"
        feature_output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/features/reprod/${run_bucket}/${run_group}/${RUN_NUMBER}"
        tt_reco_filepath="${XRD_URL}${XRD_BASEPATH}/juno/user/j/jpandre_1/tt_data_auto/${run_bucket}/${run_group}/${RUN_NUMBER}/RUN.${RUN_NUMBER}.*.EDM.user.root"
    else
        log ERROR "Unrecognized path format: ${input_file}"
        return 1  # Better to return 1 here and let caller decide to exit
    fi

    mkdir -p "${output_path}"
    mkdir -p "${reco_output_path}"
    mkdir -p "${feature_output_path}"

    log INFO "Output path: ${output_path}"
    log DEBUG "TT reco file path: ${tt_reco_filepath}"
}

#==============================
# RTRAW to ReProd filename
#==============================

rtraw_to_reprod_filename() {
    local fpath="$1"
    local fname=$(basename "${fpath}")

    local bucket_val=$(( (10#$RUN_NUMBER / 1000) * 1000 ))
    local group_val=$(( (10#$RUN_NUMBER / 100) * 100 ))
    local run_bucket=$(printf "%08d" "${bucket_val}")
    local run_group=$(printf "%08d" "${group_val}")

    local stream
    if [[ "${fname}" =~ \.([^.]+)\.[0-9]{14}\.[0-9]+_ ]]; then
        stream="${BASH_REMATCH[1]}"
    else
        stream="unknown_stream"
    fi

    local output_reprod_filename="${fname/.rtraw/.esd}"

    local base_dir="${XRD_BASEPATH}/juno/juno-reprod/${CAMPAIGN}/${stream}/${run_bucket}"

    local candidate_groups
    candidate_groups=$(xrdfs "${XRD_URL}" ls "${base_dir}" 2>/dev/null | grep -E "/${run_group}(_v[0-9]+)?/?$" | sort)

    if [[ -z "${candidate_groups}" ]]; then
        log ERROR "No run_group directory found under ${base_dir}" >&2
        exit 1
    fi

    local selected_group
    selected_group=$(basename "$(echo "${candidate_groups}" | tail -n 1)")

    local reprod_path="${XRD_URL}${XRD_BASEPATH}/juno/juno-reprod/${CAMPAIGN}/${stream}/${run_bucket}/${selected_group}/${RUN_NUMBER}/${output_reprod_filename}"

    echo "${reprod_path}"
}

#==============================
# Main Execution
#==============================

main() {
    parse_args "$@"
    load_file_lists

    source ${SOURCE_JUNOSW_PATH}

    input_file="${RTRAW_LIST[$PROC_ID]}"
    input_filename=$(basename "${input_file}")

    target_num=$(get_file_number "${input_filename}")
    if [[ -z $target_num ]]; then
        log ERROR "Cannot extract file number from ${input_filename}"
        exit 1
    fi
    target_num=$((10#$target_num))

    prev_idx=$((PROC_ID - 1))
    next_idx=$((PROC_ID + 1))

    indices_to_process=("${PROC_ID}")
    include_neighbor "${PROC_ID}" "prev"
    include_neighbor "${PROC_ID}" "next"

    indices_to_process=($(printf "%s\n" "${indices_to_process[@]}" | sort -n))
    log INFO "Files to process: ${indices_to_process[*]}"

    input_files=()
    input_correlation_files=()

    for idx in "${indices_to_process[@]}"; do
        input_corr_file="${RTRAW_LIST[$idx]}"
        local_corr_file="${TEMPDIR}/$(basename "${input_corr_file}")"
        if [[ "${USE_CORRELATION}" == "true" ]]; then
            input_main_file="${ESD_LIST[$idx]}"
        else
            input_main_file=$(rtraw_to_reprod_filename "${input_corr_file}")
        fi
        local_main_file="${TEMPDIR}/$(basename "${input_main_file}")"

        log INFO "Copying input file index ${idx}: $(basename "${local_main_file}")"
        
        xrdcp "${input_main_file}" "${local_main_file}"
        input_files+=("${local_main_file}")

        if [[ "${USE_CORRELATION}" == "true" ]]; then
            xrdcp "${input_corr_file}" "${local_corr_file}"
            input_correlation_files+=("${local_corr_file}")
        fi
    done

    resolve_output_paths "${input_file}"

    local_output_file="${TEMPDIR}/${input_filename/.rtraw/.output.root}"
    output_file="${output_path}/$(basename "${local_output_file}")"

    local_reco_output_file="${TEMPDIR}/${input_filename/.rtraw/.reco.output.root}"
    reco_output_file="${reco_output_path}/$(basename "${local_reco_output_file}")"

    local_feature_output_file="${TEMPDIR}/${input_filename/.rtraw/.feature.output.root}"
    feature_output_file="${feature_output_path}/$(basename "${local_feature_output_file}")"


    if (( prev_idx >= 0 )); then
        if [[ "${USE_CORRELATION}" == "true" ]]; then
            prev_file_local="$(basename "${ESD_LIST[${prev_idx}]}")"
        else
            prev_file_local="$(basename "$(rtraw_to_reprod_filename "${RTRAW_LIST[${prev_idx}]}")")"
        fi
    else
        prev_file_local=""
    fi

    if (( next_idx < ${#RTRAW_LIST[@]} )); then
        if [[ "${USE_CORRELATION}" == "true" ]]; then
            next_file_local="$(basename "${ESD_LIST[${next_idx}]}")"
        else
            next_file_local="$(basename "$(rtraw_to_reprod_filename "${RTRAW_LIST[${next_idx}]}")")"
        fi
    else
        next_file_local=""
    fi

    log INFO "Context previous file: ${prev_file_local:-<none>}"
    log INFO "Context next file: ${next_file_local:-<none>}"

    log INFO "Running run.py..."

    python_args=(
        "run.py"
        "--input" "${input_files[@]}"
    )

    if [[ "${USE_CORRELATION}" == "true" ]]; then
        python_args+=("--input-correlation" "${input_correlation_files[@]}")
    fi

    python_args+=(
        "--output" "${local_output_file}"
        "--context-previous-filename" "${prev_file_local}"
        "--context-next-filename" "${next_file_local}"
        "--tt-reco-filepath" "${tt_reco_filepath}"
        "--reco-output" "${local_reco_output_file}"
        "--feature-output" "${local_feature_output_file}"
        "--cluster" "${CLUSTER}"
        "${EXTRA_ARGS[@]}"
    )

    python "${python_args[@]}"

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

    if cp "${local_feature_output_file}" "${feature_output_file}"; then
        log INFO "Output copied to ${feature_output_file}"
    else
        log ERROR "Failed to copy ${local_feature_output_file} to ${feature_output_file}"
        exit 1
    fi

    log INFO "Completed run ${RUN_NUMBER} (PROC_ID=${PROC_ID})"
}

main "$@"