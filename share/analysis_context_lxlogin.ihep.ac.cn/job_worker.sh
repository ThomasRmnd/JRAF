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

RTRAW_LIST_FILE="$LIST_BASE/rtraw_list/run_${RUN_NUMBER}.txt"
ESD_LIST_FILE="$LIST_BASE/esd_list/run_${RUN_NUMBER}.txt"

mapfile -t rtraw_list < <(xrdfs "$EOS_BASE" cat "$RTRAW_LIST_FILE")
mapfile -t esd_list < <(xrdfs "$EOS_BASE" cat "$ESD_LIST_FILE")

if (( ${#rtraw_list[@]} != ${#esd_list[@]} )); then
    log "Error: mismatch in rtraw (${#rtraw_list[@]}) vs esd (${#esd_list[@]})"
    exit 1
fi

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
input_track_path=""
if [[ "$input_esd_file" =~ /eos/juno/esd/([^/]+)/([0-9]{4})/([0-9]{4})/RUN\.([0-9]+)\. ]]; then
    esd_version="${BASH_REMATCH[1]}"
    year="${BASH_REMATCH[2]}"
    monthday="${BASH_REMATCH[3]}"
    output_path="/junofs/users/traymond/analysis/ibd/${year}/${monthday}"
    input_track_path="/scratchfs/juno/traymond/analysis/ibd/${year}/${monthday}"
elif [[ "$input_esd_file" =~ /eos/juno-kup/([^/]+)/([^/]+)/([^/]+)/([^/]+)/([^/]+)/RUN\.([0-9]+)\. ]]; then
    campaign="${BASH_REMATCH[1]}"
    stream="${BASH_REMATCH[2]}"
    run_bucket="${BASH_REMATCH[3]}"
    run_group="${BASH_REMATCH[4]}"
    run_number="${BASH_REMATCH[5]}"
    output_path="/junofs/users/traymond/analysis/ibd/${run_bucket}/${run_group}/${RUN_NUMBER}"
    input_track_path="/scratchfs/juno/traymond/analysis/ibd/${run_bucket}/${run_group}/${RUN_NUMBER}"
else
    echo "Error: unrecognized esd file path format: $input_esd_file"
    exit 1
fi
mkdir -p "$output_path"

echo "output_path: ${output_path}" 1>&2

input_track_filename="${input_esd_filename/.esd/.track.rec}"
input_track_file="$input_track_path/$input_track_filename"

target_num=$(get_file_number "$input_esd_filename")

local_output_file="$TEMP/${input_esd_filename/.esd/.output.root}"
output_file="$output_path/$(basename "$local_output_file")"

if [[ -z $target_num ]]; then
    log "Error: cannot extract file number from $input_track_filename"
    exit 1
fi
target_num=$((10#$target_num))

log "Target file: $input_track_filename (index=$PROC_ID, number=$target_num)"

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

track_files=()
rtraw_files=()
idx=0

for i in "${indices_to_process[@]}"; do
    esd_file="${esd_list[$i]}"
    rtraw_file="${rtraw_list[$i]}"
    fname=$(basename "$esd_file")
    log "Processing file: $fname"

    local_esd="$TEMP/$fname"
    local_rtraw="$TEMP/$(basename "$rtraw_file")"
    track_file="$input_track_path/${fname/.esd/.track.rec}"
    local_track="$TEMP/${fname/.esd/.track.rec}"

    xrdcp "$esd_file" "$local_esd"
    xrdcp "$rtraw_file" "$local_rtraw"
    cp "$track_file" "$local_track"

    track_files+=("$local_track")
    rtraw_files+=("$local_rtraw")
    rtraw_files+=("$local_esd")
    idx=$((idx + 1))
done

if (( idx == 0 )); then
    log "ERROR: No files processed"
    exit 1
fi

log "Running run_analysis.py on ${#track_files[@]} files"
python run_analysis.py \
    --input "${track_files[@]}" \
    --input-rtraw "${rtraw_files[@]}" \
    --output "$local_output_file" \
    --target-input-filename "$input_track_filename" \
    "${EXTRA_ARGS[@]}"

cp "$local_output_file" "$output_file"
log "Final result copied to $output_file"
log "Done for run $RUN_NUMBER (target $input_track_filename)"