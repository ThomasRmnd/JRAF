#!/bin/bash
set -e

PROC_ID="$1"
RUN_NUMBER="$2"
LIST_BASE="$3"
OUTPUT_DIR="$4"

shift 4

EOS_BASE="root://junoeos01.ihep.ac.cn/"
SKIP_IF_EXISTS=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-if-exists)
            SKIP_IF_EXISTS=true
            shift
            ;;
        *) 
            EXTRA_ARGS+=("$1")
            shift
            ;;
    esac
done

RTRAW_LIST_FILE="$LIST_BASE/rtraw_list/run_${RUN_NUMBER}.txt"
ESD_LIST_FILE="$LIST_BASE/esd_list/run_${RUN_NUMBER}.txt"

input_rtraw_file=$(xrdfs "$EOS_BASE" cat "$RTRAW_LIST_FILE" | sed -n "$((PROC_ID + 1))p")
input_esd_file=$(xrdfs "$EOS_BASE" cat "$ESD_LIST_FILE" | sed -n "$((PROC_ID + 1))p")

if [[ -z "$input_rtraw_file" || -z $input_esd_file ]]; then
    echo "No input file found for PROC_ID=$PROC_ID"
    exit 1
fi

input_rtraw_filename=$(basename "$input_rtraw_file")
input_esd_filename=$(basename "$input_esd_file")
local_input_rtraw_file="$TEMP/$input_rtraw_filename"
local_input_esd_file="$TEMP/$input_esd_filename"
local_track_file="$TEMP/${input_esd_filename/.esd/.track.rec}"
local_output_file="$TEMP/${input_esd_filename/.esd/.output.root}"
output_file="$OUTPUT_DIR/$(basename "$local_output_file")"

echo "Output filename: $output_file"

if [[ "$SKIP_IF_EXISTS" == true && -f "$output_file" ]]; then
    echo "Skipping PROC_ID=$PROC_ID — output already exists: $output_file"
    exit 0
fi

echo "Copying rtraw file from EOS: $local_input_rtraw_file"
xrdcp "${input_rtraw_file}" "$local_input_rtraw_file"
echo "Copying esd file from EOS: $local_input_esd_file"
xrdcp "${input_esd_file}" "$local_input_esd_file"

source /afs/ihep.ac.cn/users/t/traymond/J25.3.0/git_junosw_J25_load.sh
echo "TUTORIALROOT = ${TUTORIALROOT}"

echo "Running run_muon_wrtraw.py for file: $local_input_esd_file"

python run_muon_wrtraw.py \
  --input "$local_input_esd_file" \
  --input-rtraw "$local_input_rtraw_file" \
  --output "$local_track_file" \
  "${EXTRA_ARGS[@]}"

echo "Running run_analysis_wrtraw.py for file: $local_track_file"

python run_analysis_wrtraw.py \
  --input "$local_track_file" \
  --input-rtraw "$local_input_rtraw_file" \
  --output "$local_output_file" \
  "${EXTRA_ARGS[@]}"

echo "Copying result to $output_file"
cp "$local_output_file" "$output_file"

echo "Done (PROC_ID=$PROC_ID)"