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

source /pbs/home/t/traymond/share/bash/logging.sh
export X509_USER_PROXY=/sps/juno/jdeandre/rtraw_ThomasRaymond/.cert_traymond_juno_user

#==============================
# Configuration defaults
#==============================

XRD_URL="root://xrootd-archive.cr.cnaf.infn.it:1095/"

#==============================
# Parse command-line arguments
#==============================

usage() {
    cat <<EOF
Usage: $(basename "$0") <run_number>

Compare reconstruction with TT reconstruction for a ReProd summary file identified by run and process ID.

Arguments:
  <run_number>    Run number to process
EOF
}

parse_args() {
    if (( $# < 1 )); then
        log ERROR "Missing required arguments"
        usage >&2
        exit 1
    fi

    RUN_NUMBER="$1"
}

#==============================
# Main Execution
#==============================

main() {
    parse_args "$@"

    input_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/summary"
    input_filename="RUN.${RUN_NUMBER}.output.reprod25c.cca.root"
    input_file="${input_path}/${input_filename}"
    local_input_file="${TMPDIR}/${input_filename}"

    bucket_val=$(( (10#$RUN_NUMBER / 1000) * 1000 ))
    group_val=$(( (10#$RUN_NUMBER / 100) * 100 ))
    run_bucket=$(printf "%08d" "$bucket_val")
    run_group=$(printf "%08d" "$group_val")
    input_tt_path="${XRD_URL}/production/storm/dirac/juno/juno-reprod/TT25A/J25.4.3-patched/user_rec/${run_bucket}/${run_group}/${RUN_NUMBER}/RUN.${RUN_NUMBER}.*"

    output_path="/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/comparison"
    output_filename="RUN.${RUN_NUMBER}.output.reprod25c.cca.root"
    output_file="${output_path}/${output_filename}"
    local_output_file="${TMPDIR}/${output_filename}"

    log INFO "Copying input from ${input_file} to ${local_input_file}"
    cp "${input_file}" "${local_input_file}"

    log INFO "Running compare_with_TT.C..."
    root -l -b -q "compare_with_TT.C(\"$local_input_file\", \"$input_tt_path\", \"$local_output_file\")"

    mkdir -p "${output_path}"
    log INFO "Copying output from ${local_output_file} to ${output_file}"
    cp "${local_output_file}" "${output_file}"

    log INFO "Completed run ${RUN_NUMBER}"
}

main "$@"