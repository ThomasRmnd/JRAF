#!/bin/bash

#--------------------------------------------------------------------------------------------------
#  JUNO Job Worker
#  Purpose: Process ReProd
#--------------------------------------------------------------------------------------------------

set -euo pipefail
IFS=$'\n\t'

#==============================
# Utility functions
#==============================

HOSTNAME=$(hostname -f 2>/dev/null || hostname)
if [[ "${HOSTNAME}" =~ ^cca[0-9]+\.in2p3\.fr$ ]]; then
    # Detect CC-IN2P3 cluster
    CLUSTER="CC-IN2P3"
    source /pbs/home/t/traymond/share/bash/logging.sh
    TEMPDIR=${TMPDIR}
    SOURCE_JUNOSW_PATH="/pbs/home/t/traymond/J25.6.1_Modified/git_junosw_J25_load.sh"
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
Usage: $(basename "$0") <run_number> <site>

Process a single file identified by run and process ID.

Arguments:
  <run_number>    Run number to process
  <site>          Storage site selection {EOS|CNAF}
EOF
}

parse_args() {
    if (( $# < 2 )); then
        echo "ERROR: Missing arguments for CC-IN2P3" >&2
        usage >&2
        exit 1
    fi

    RUN_NUMBER="$1"; shift
    SITE="$1"; shift

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
}

#==============================
# Path Handling
#==============================

resolve_input_paths() {
    local rtraw_list_file="${LIST_BASE}/rtraw_list/run_${RUN_NUMBER}.txt"
    log INFO "Listing ROOT files from EOS..."
    mapfile -t RTRAW_LIST < <(xrdfs "${XRD_URL_EOS}" cat "${rtraw_list_file}")
    log INFO "Number of RTRAW files: ${#RTRAW_LIST[@]}"
    input_file="${RTRAW_LIST[0]}"

    if [[ "${input_file}" =~ /juno/esd/([^/]+)/([0-9]{4})/([0-9]{4})/RUN\.${RUN_NUMBER}\. ]]; then
        local esd_version="${BASH_REMATCH[1]}"
        local year="${BASH_REMATCH[2]}"
        local monthday="${BASH_REMATCH[3]}"

        tt_reco_filepath="${XRD_URL}${XRD_BASEPATH}/juno/user/j/jpandre_1/tt_data_auto/${year}/${monthday}/RUN.${RUN_NUMBER}.*.EDM.user.root"
    elif [[ "${input_file}" =~ /juno/juno-([^/]+)/([^/]+)/([^/]+)/([^/]+)/([^/]+)/${RUN_NUMBER}/RUN\.${RUN_NUMBER}\. ]]; then
        local type="${BASH_REMATCH[1]}"
        local campaign="${BASH_REMATCH[2]}"
        local stream="${BASH_REMATCH[3]}"
        local run_bucket="${BASH_REMATCH[4]}"
        local run_group="${BASH_REMATCH[5]}"

        elif (( RUN_NUMBER >= 10176 && RUN_NUMBER <= 10479 )); then
            tt_reco_filepath="${XRD_URL}${XRD_BASEPATH}/juno/user/j/jpandre_1/tt_data_auto/${run_bucket}/${run_group}/${RUN_NUMBER}/RUN.${RUN_NUMBER}.*.EDM.user.root"
        elif (( RUN_NUMBER >= 10480 )); then
            tt_reco_filepath="${XRD_URL}${XRD_BASEPATH}/juno/juno-reprod/TT25A/J25.4.3-patched/user_rec/${run_bucket}/${run_group}/${RUN_NUMBER}/RUN.${RUN_NUMBER}.*.EDM.user.root"
        else
            log ERROR "No TT reco path rule defined for run ${RUN_NUMBER}"
            exit 1
        fi
    else
        log ERROR "Unrecognized path format: ${input_file}"
        return 1  # Better to return 1 here and let caller decide to exit
    fi

    log INFO "TT-Reco filepath resolved: ${tt_reco_filepath}"
}

#==============================
# Main Execution
#==============================

main() {
    parse_args "$@"
    resolve_input_paths

    input_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/summary"
    input_filename="RUN.${RUN_NUMBER}.output.reprod25c.cca.root"
    input_file="${input_path}/${input_filename}"
    local_input_file="${TMPDIR}/${input_filename}"

    output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/comparison"
    output_filename="RUN.${RUN_NUMBER}.output.reprod25c.cca.root"
    output_file="${output_path}/${output_filename}"
    local_output_file="${TMPDIR}/${output_filename}"

    log INFO "Copying input from ${input_file} to ${local_input_file}"
    cp "${input_file}" "${local_input_file}"

    log INFO "Running compare_with_TT.C..."
    root -l -b -q "compare_with_TT.C(\"${local_input_file}\", \"${tt_reco_filepath}\", \"${local_output_file}\")"

    mkdir -p "${output_path}"
    log INFO "Copying output from ${local_output_file} to ${output_file}"
    cp "${local_output_file}" "${output_file}"

    log INFO "Completed run ${RUN_NUMBER}"
}

main "$@"