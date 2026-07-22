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

INPUT_SUFFIX="reconstruction.root"
OUTPUT_SUFFIX="reconstruction.cdwpttchi2.root"

#==============================
# Parse command-line arguments
#==============================

usage() {
    cat <<EOF
Usage: $(basename "$0") <run_number>

Process a single file identified by run and process ID.

Arguments:
  <run_number>    Run number to process
EOF
}

parse_args() {
    if (( $# < 2 )); then
        echo "ERROR: Missing arguments for CC-IN2P3" >&2
        usage >&2
        exit 1
    fi

    RUN_NUMBER="$1"; shift
}

#==============================
# Main Execution
#==============================

main() {
    parse_args "$@"

    source ${SOURCE_JUNOSW_PATH}

    input_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/summary"
    input_filename="RUN.${RUN_NUMBER}.${INPUT_SUFFIX}"
    input_file="${input_path}/${input_filename}"
    local_input_file="${TEMPDIR}/${input_filename}"

    output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/summary/CdWpTtChi2"
    output_filename="RUN.${RUN_NUMBER}.${OUTPUT_SUFFIX}"
    output_file="${output_path}/${output_filename}"
    local_output_file="${TEMPDIR}/${output_filename}"

    log INFO "Copying input from ${input_file} to ${local_input_file}"
    cp "${input_file}" "${local_input_file}"

    log INFO "Running extract_cdwpttchi2.C..."
    root -l -b -q "extract_cdwpttchi2.C(\"${local_input_file}\", \"${local_output_file}\")"

    mkdir -p "${output_path}"
    log INFO "Copying output from ${local_output_file} to ${output_file}"
    cp "${local_output_file}" "${output_file}"

    log INFO "Completed run ${RUN_NUMBER}"
}

main "$@"