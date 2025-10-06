#!/bin/bash
set -e

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

PROC_ID="$1"
RUN_NUMBER="$2"
LIST_BASE="$3"
OUTPUT_DIR="$4"

shift 4
EXTRA_ARGS=("$@")

EOS_BASE="root://junoeos01.ihep.ac.cn/"

RTRAW_LIST_FILE="$LIST_BASE/rtraw_list/run_${RUN_NUMBER}.txt"
ESD_LIST_FILE="$LIST_BASE/esd_list/run_${RUN_NUMBER}.txt"

mapfile -t rtraw_list < <(xrdfs "$EOS_BASE" cat "$RTRAW_LIST_FILE")
mapfile -t esd_list   < <(xrdfs "$EOS_BASE" cat "$ESD_LIST_FILE")

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

target_esd="${esd_list[$PROC_ID]}"
target_fname=$(basename "$target_esd")
target_num=$(get_file_number "$target_fname")
if [[ -z $target_num ]]; then
    log "Error: cannot extract file number from $target_fname"
    exit 1
fi
target_num=$((10#$target_num))

log "Target file: $target_fname (index=$PROC_ID, number=$target_num)"

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

trackfile_path="${OUTPUT_DIR/\/junofs\/users/\/scratchfs\/juno}"

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
    track_file="$trackfile_path/${fname/.esd/.track.rec}"
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

first_input=$(basename "${esd_list[${indices_to_process[0]}]}")
prefix=$(echo "$first_input" | sed -E 's/^(.*\.)[0-9]{14}\..*/\1/')
suffix=$(echo "$first_input" | sed -E 's/^.*\.[0-9]{14}\.[0-9]+(.*)$/\1/')
final_output="$TEMP/${prefix}${target_num}${suffix}"

log "Running run_analysis.py on ${#track_files[@]} files"
python run_analysis.py \
    --input "${track_files[@]}" \
    --input-rtraw "${rtraw_files[@]}" \
    --output "$final_output" \
    --target-input-filename "$target_fname" \
    "${EXTRA_ARGS[@]}"

cp "$final_output" "$OUTPUT_DIR/"
log "Final result copied to $OUTPUT_DIR/$(basename "$final_output")"
log "Done for run $RUN_NUMBER (target $target_fname)"