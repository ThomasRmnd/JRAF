#!/bin/bash
set -e

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

prefix=$(echo "$first_input" | sed -E 's/^(.*\.[0-9]{14})\..*/\1/')
suffix=$(echo "$first_input" | sed -E 's/^.*\.[0-9]{14}\.[0-9]+(.*)$/\1/')

final_output="$TEMP/${prefix}.${RANGE_START}-${RANGE_END}${suffix}"
echo "Output file: $final_output"

track_files=()
rtraw_files=()

first_flag="--first-reconstruction-file True"
idx=0

source /afs/ihep.ac.cn/users/t/traymond/J25.3.0/git_junosw_J25_load.sh
echo "TUTORIALROOT = ${TUTORIALROOT}"

for ((i=0; i<${#esd_list[@]}; i++)); do
    esd_file="${esd_list[$i]}"
    rtraw_file="${rtraw_list[$i]}"

    fname=$(basename "$esd_file")
    file_number=$(echo "$fname" | sed -n 's/^.*\.[0-9]\{14\}\.\([0-9]*\)_.*$/\1/p')
    file_number=$((file_number + 0))

    if (( file_number < RANGE_START || file_number > RANGE_END )); then
        continue
    fi

    # Local names
    local_esd="$TEMP/$(basename "$esd_file")"
    local_rtraw="$TEMP/$(basename "$rtraw_file")"
    local_track="$TEMP/${fname/.esd/.track.rec}"

    echo "Copying files for file_number=$file_number"
    xrdcp "$esd_file" "$local_esd"
    xrdcp "$rtraw_file" "$local_rtraw"

    echo "Running run_muon.py on $fname"
    if (( idx == 0 )); then
        python run_muon.py \
          --input "$local_esd" \
          --input-rtraw "$local_rtraw" \
          --output "$local_track" \
          --first-reconstruction-file True \
          "${EXTRA_ARGS[@]}"
    else
        python run_muon.py \
          --input "$local_esd" \
          --input-rtraw "$local_rtraw" \
          --output "$local_track" \
          --first-reconstruction-file False \
          "${EXTRA_ARGS[@]}"
    fi

    track_files+=("$local_track")
    rtraw_files+=("$local_rtraw")
    ((idx++))
done

if (( idx == 0 )); then
    echo "No files found in range $RANGE_START-$RANGE_END"
    exit 1
fi

echo "Running run_analysis.py on ${#track_files[@]} files"

python run_analysis.py \
  --input "${track_files[@]}" \
  --input-rtraw "${rtraw_files[@]}" \
  --output "$final_output" \
  "${EXTRA_ARGS[@]}"

cp "$final_output" "$OUTPUT_DIR/"
echo "Final result copied to $OUTPUT_DIR/$(basename "$final_output")"

echo "Done for run $RUN_NUMBER, range $RANGE_START-$RANGE_END"