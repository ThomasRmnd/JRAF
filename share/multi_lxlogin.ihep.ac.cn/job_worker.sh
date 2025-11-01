#!/bin/bash
set -e

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

RANGE_START="$1"
RANGE_END="$2"
RUN_NUMBER="$3"
LIST_BASE="$4"
shift 4
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

first_input="${esd_list[0]}"
first_input_filename=$(basename "$first_input")

output_path=""
trackfile_path=""
if [[ "$first_input" =~ /eos/juno/esd/([^/]+)/([0-9]{4})/([0-9]{4})/RUN\.([0-9]+)\. ]]; then
    esd_version="${BASH_REMATCH[1]}"
    year="${BASH_REMATCH[2]}"
    monthday="${BASH_REMATCH[3]}"
    output_path="/junofs/users/traymond/analysis/ibd/${year}/${monthday}"
    trackfile_path="/scratchfs/juno/traymond/analysis/ibd/${year}/${monthday}"
elif [[ "$first_input" =~ /eos/juno-kup/([^/]+)/([^/]+)/([^/]+)/([^/]+)/([^/]+)/RUN\.([0-9]+)\. ]]; then
    campaign="${BASH_REMATCH[1]}"
    stream="${BASH_REMATCH[2]}"
    run_bucket="${BASH_REMATCH[3]}"
    run_group="${BASH_REMATCH[4]}"
    run_number="${BASH_REMATCH[5]}"
    output_path="/junofs/users/traymond/analysis/ibd/${run_bucket}/${run_group}/${RUN_NUMBER}"
    trackfile_path="/scratchfs/juno/traymond/analysis/ibd/${run_bucket}/${run_group}/${RUN_NUMBER}"
else
    echo "Error: unrecognized esd file path format: $first_input"
    exit 1
fi
mkdir -p "$output_path"

echo "output_path: ${output_path}" 1>&2
echo "trackfile_path: ${trackfile_path}" 1>&2

prefix=$(echo "$first_input_filename" | sed -E 's/^(.*\.)[0-9]{14}\..*/\1/')
suffix=$(echo "$first_input_filename" | sed -E 's/^.*\.[0-9]{14}\.[0-9]+(.*)$/\1/')
date="${timestamp:0:8}"
final_output="$TEMP/${prefix}${date}.${RANGE_START}-${RANGE_END}${suffix}"
echo "Output file: $final_output"

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
    
    log "Copying ESD file $esd_file"
    xrdcp "$esd_file" "$local_esd"
    log "ESD $local_esd copied"
    
    log "Copying RTRAW file $rtraw_file"
    xrdcp "$rtraw_file" "$local_rtraw"
    log "RTRAW $local_rtraw copied"

    log "Copying track file $track_file"
    cp "$track_file" "$local_track"
    log "Track $local_track copied"
    
    track_files+=("$local_track")
    rtraw_files+=("$local_rtraw")
    rtraw_files+=("$local_esd")
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

cp "$final_output" "$output_path/"
log "Final result copied to $output_path/$(basename "$final_output")"
log "Done for run $RUN_NUMBER, range $RANGE_START-$RANGE_END"