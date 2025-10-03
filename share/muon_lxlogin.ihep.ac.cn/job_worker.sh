#!/bin/bash
set -e

PROC_ID="$1"
RANGE_START="$2"
RANGE_END="$3"
RUN_NUMBER="$4"
LIST_BASE="$5"
OUTPUT_DIR="$6"
shift 6
EXTRA_ARGS=("$@")

EOS_BASE="root://junoeos01.ihep.ac.cn/"
RTRAW_LIST_FILE="$LIST_BASE/rtraw_list/run_${RUN_NUMBER}.txt"
ESD_LIST_FILE="$LIST_BASE/esd_list/run_${RUN_NUMBER}.txt"

FILE_IDX=$((RANGE_START + PROC_ID))

if (( FILE_NUMBER > RANGE_END )); then
    echo "ERROR: FILE_NUMBER=$FILE_NUMBER exceeds RANGE_END=$RANGE_END"
    exit 1
fi

input_rtraw_file=$(xrdfs "$EOS_BASE" cat "$RTRAW_LIST_FILE" | sed -n "$((FILE_IDX))p")
input_esd_file=$(xrdfs "$EOS_BASE" cat "$ESD_LIST_FILE" | sed -n "$((FILE_IDX))p")

if [[ -z "$input_rtraw_file" || -z "$input_esd_file" ]]; then
    echo "No input file found for FILE_IDX=$FILE_IDX, PROC_ID=$PROC_ID"
    exit 1
fi

input_rtraw_filename=$(basename "$input_rtraw_file")
input_esd_filename=$(basename "$input_esd_file")

local_input_rtraw_file="$TEMP/$input_rtraw_filename"
local_input_esd_file="$TEMP/$input_esd_filename"
local_output_file="$TEMP/${input_esd_filename/.esd/.track.rec}"
output_file="$OUTPUT_DIR/$(basename "$local_output_file")"

echo "Output filename: $output_file"

echo "Copying rtraw file from EOS: $local_input_rtraw_file"
xrdcp "${input_rtraw_file}" "$local_input_rtraw_file"

echo "Copying esd file from EOS: $local_input_esd_file"
xrdcp "${input_esd_file}" "$local_input_esd_file"

source /afs/ihep.ac.cn/users/t/traymond/J25.3.0/git_junosw_J25_load.sh
echo "TUTORIALROOT = ${TUTORIALROOT}"

echo "Running run_muon.py for file: $local_input_esd_file (FILE_IDX=$FILE_IDX, PROC_ID=$PROC_ID)"

if (( PROC_ID == 0 )); then
    echo "Using --first-reconstruction-file flag"
    python run_muon.py \
        --input "$local_input_esd_file" \
        --input-rtraw "$local_input_rtraw_file" \
        --output "$local_output_file" \
        --first-reconstruction-file \
        "${EXTRA_ARGS[@]}"
else
    echo "Not using --first-reconstruction-file flag"
    python run_muon.py \
        --input "$local_input_esd_file" \
        --input-rtraw "$local_input_rtraw_file" \
        --output "$local_output_file" \
        "${EXTRA_ARGS[@]}"
fi

echo "Copying result to $output_file"
cp "$local_output_file" "$output_file"

echo "Done (FILE_IDX=$FILE_IDX, PROC_ID=$PROC_ID)"