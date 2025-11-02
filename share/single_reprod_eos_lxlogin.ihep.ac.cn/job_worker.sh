#!/bin/bash
set -e

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

PROC_ID="$1"
RUN_NUMBER="$2"
LIST_BASE="$3"
shift 3
EXTRA_ARGS=("$@")

EOS_BASE="root://junoeos01.ihep.ac.cn/"

ESD_LIST_FILE="$LIST_BASE/esd_list/run_${RUN_NUMBER}.txt"

mapfile -t esd_list < <(xrdfs "$EOS_BASE" cat "$ESD_LIST_FILE")

num_files=${#esd_list[@]}

get_file_number() {
    local fname=$1
    if [[ $fname =~ \.[0-9]{14}\.([0-9]+)_ ]]; then
        echo "${BASH_REMATCH[1]}"
    else
        echo ""
    fi
}

input_esd_file="${esd_list[$PROC_ID]}"
input_esd_filename=$(basename "$input_esd_file")

output_path=""
tt_reco_filepath=""
if [[ "$input_esd_file" =~ /eos/juno/esd/([^/]+)/([0-9]{4})/([0-9]{4})/RUN\.([0-9]+)\. ]]; then
    esd_version="${BASH_REMATCH[1]}"
    year="${BASH_REMATCH[2]}"
    monthday="${BASH_REMATCH[3]}"
    output_path="/junofs/users/traymond/analysis/ibd/${year}/${monthday}"
    tt_reco_filepath="${EOS_BASE}/eos/juno/dirac/juno/user/j/jpandre_1/tt_data_auto/${year}/${monthday}/RUN.${RUN_NUMBER}.*.EDM.user.root"
elif [[ "$input_esd_file" =~ /eos/juno-kup/([^/]+)/([^/]+)/([^/]+)/([^/]+)/([^/]+)/RUN\.([0-9]+)\. ]]; then
    campaign="${BASH_REMATCH[1]}"
    stream="${BASH_REMATCH[2]}"
    run_bucket="${BASH_REMATCH[3]}"
    run_group="${BASH_REMATCH[4]}"
    run_number="${BASH_REMATCH[5]}"
    output_path="/junofs/users/traymond/analysis/ibd/${run_bucket}/${run_group}/${RUN_NUMBER}"
    tt_reco_filepath="${EOS_BASE}/eos/juno/dirac/juno/user/j/jpandre_1/tt_data_auto/${run_bucket}/${run_group}/${RUN_NUMBER}/RUN.${RUN_NUMBER}.*.EDM.user.root"
else
    log "Error: unrecognized esd file path format: $input_esd_file"
    exit 1
fi
mkdir -p "$output_path"

log "output_path: ${output_path}" 1>&2
log "tt_reco_filepath: $tt_reco_filepath" 1>&2

local_input_esd_file="$TEMP/$input_esd_filename"
local_output_file="$TEMP/${input_esd_filename/.esd/.output.root}"
output_file="$output_path/$(basename "$local_output_file")"

target_num=$(get_file_number "$input_esd_filename")

if [[ -z $target_num ]]; then
    log "Error: cannot extract file number from $input_esd_filename"
    exit 1
fi
target_num=$((10#$target_num))

indices_to_process=("$PROC_ID")

# Previous file (if exists and number = target_num - 1)
if (( PROC_ID > 0 )); then
    prev_fname=$(basename "${esd_list[$((PROC_ID - 1))]}")
    prev_num=$(get_file_number "$prev_fname")
    prev_num=$((10#$prev_num))
    if (( prev_num == target_num - 1 )); then
        indices_to_process=("$((PROC_ID - 1))" "${indices_to_process[@]}")
        log "Including previous file: $prev_fname (number=$prev_num)"
    fi
fi

# Next file (if exists and number = target_num + 1)
if (( PROC_ID < num_files - 1 )); then
    next_fname=$(basename "${esd_list[$((PROC_ID + 1))]}")
    next_num=$(get_file_number "$next_fname")
    next_num=$((10#$next_num))
    if (( next_num == target_num + 1 )); then
        indices_to_process+=("$((PROC_ID + 1))")
        log "Including next file: $next_fname (number=$next_num)"
    fi
fi

log "Files to process (indices): ${indices_to_process[*]}"

# --- Environment setup ---
source /afs/ihep.ac.cn/users/t/traymond/J25.3.0/git_junosw_J25_load.sh
log "TUTORIALROOT = ${TUTORIALROOT}"

esd_files=()
idx=0

for i in "${indices_to_process[@]}"; do
    esd_file="${esd_list[$i]}"
    fname=$(basename "$esd_file")
    log "Processing file: $fname"

    local_esd="$TEMP/$fname"

    xrdcp "$esd_file" "$local_esd"

    esd_files+=("$local_esd")
    idx=$((idx + 1))
done

if (( idx == 0 )); then
    log "ERROR: No files processed"
    exit 1
fi

log "Running run.py on ${#esd_files[@]} files"
python run.py \
    --input "${esd_files[@]}" \
    --output "$local_output_file" \
    --target-input-filename "$input_esd_filename" \
    --tt-reco-filepath "$tt_reco_filepath" \
    "${EXTRA_ARGS[@]}"

log "Copying result to $output_file"
cp "$local_output_file" "$output_file"

log "Done for run ${RUN_NUMBER} (PROC_ID=${PROC_ID})"