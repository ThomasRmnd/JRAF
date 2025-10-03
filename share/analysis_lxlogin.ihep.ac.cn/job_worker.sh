#!/bin/bash
set -e

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

RANGE_START="$1"
RANGE_END="$2"
RUN_NUMBER="$3"
LIST_BASE="$4"
OUTPUT_DIR="$5"

shift 5
EXTRA_ARGS=("$@")

EOS_BASE="root://junoeos01.ihep.ac.cn/"

RTRAW_LIST_FILE="$LIST_BASE/rtraw_list/run_${RUN_NUMBER}.txt"
ESD_LIST_FILE="$LIST_BASE/esd_list/run_${RUN_NUMBER}.txt"

mapfile -t rtraw_list < <(xrdfs "$EOS_BASE" cat "$RTRAW_LIST_FILE")
mapfile -t esd_list   < <(xrdfs "$EOS_BASE" cat "$ESD_LIST_FILE")

if (( ${#rtraw_list[@]} != ${#esd_list[@]} )); then
    echo "Error: mismatch in rtraw (${#rtraw_list[@]}) vs esd (${#esd_list[@]})"
    exit 1
fi

echo "Processing run $RUN_NUMBER, range $RANGE_START-$RANGE_END"
echo "Temporary dir: $TEMP"

first_input=$(basename "${esd_list[0]}")

prefix=$(echo "$first_input" | sed -E 's/^(.*\.)[0-9]{14}\..*/\1/')
suffix=$(echo "$first_input" | sed -E 's/^.*\.[0-9]{14}\.[0-9]+(.*)$/\1/')
date="${timestamp:0:8}"
final_output="$TEMP/${prefix}${date}.${RANGE_START}-${RANGE_END}${suffix}"
echo "Output file: $final_output"

trackfile_path="${output_path/\/junofs\/users/\/scratchfs\/juno}"

track_files=()
rtraw_files=()
idx=0

source /afs/ihep.ac.cn/users/t/traymond/J25.3.0/git_junosw_J25_load.sh
echo "TUTORIALROOT = ${TUTORIALROOT}"

for ((i=0; i<${#esd_list[@]}; i++)); do
    log "=== Loop iteration $i of ${#esd_list[@]} ==="
    
    esd_file="${esd_list[$i]}"
    rtraw_file="${rtraw_list[$i]}"
    fname=${esd_file##*/}
    log "Current file: $fname"
    
    if [[ $fname =~ \.[0-9]{14}\.([0-9]+)_ ]]; then
        file_number=${BASH_REMATCH[1]}
        file_number=$((10#$file_number))
    else
        log "Warning: could not extract file number from $fname"
        continue
    fi
    
    log "File number: $file_number (range: $RANGE_START-$RANGE_END)"
    
    if (( file_number < RANGE_START )); then
        log "Skipping file_number=$file_number (below range)"
        continue
    fi
    
    if (( file_number > RANGE_END )); then
        log "Stopping at file_number=$file_number (above range)"
        break
    fi
    
    log "Processing file_number=$file_number"
    
    local_esd="$TEMP/$fname"
    local_rtraw="$TEMP/$(basename "$rtraw_file")"
    track_file="$trackfile_path/${fname/.esd/.track.rec}"
    local_track="$TEMP/${fname/.esd/.track.rec}"
    
    log "Copying ESD file..."
    xrdcp "$esd_file" "$local_esd"
    log "ESD copied"
    
    log "Copying RTRAW file..."
    xrdcp "$rtraw_file" "$local_rtraw"
    log "RTRAW copied"

    log "Copying track file..."
    cp "$track_file" "$local_track"
    log "Track copied"
    
    track_files+=("$local_esd")
    rtraw_files+=("$local_rtraw")
    rtraw_files+=("$local_track")
    idx=$((idx + 1))
    
    log "=== Completed iteration $i, idx now=$idx ==="
done

log "Loop finished. Total files processed: $idx"
log "Track files: ${track_files[*]}"

if (( idx == 0 )); then
    log "ERROR: No files found in range $RANGE_START-$RANGE_END"
    exit 1
fi

log "Running run_analysis.py on ${#track_files[@]} files"
python run_analysis.py \
    --input "${track_files[@]}" \
    --input-rtraw "${rtraw_files[@]}" \
    --output "$final_output" \
    "$@"
log "run_analysis.py completed"

cp "$final_output" "$OUTPUT_DIR/"
log "Final result copied to $OUTPUT_DIR/$(basename "$final_output")"
log "Done for run $RUN_NUMBER, range $RANGE_START-$RANGE_END"