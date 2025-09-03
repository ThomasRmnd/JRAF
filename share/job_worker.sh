#!/bin/bash
set -e

PROC_ID="$1"
LIST_FILE="$2"
OUTPUT_DIR="$3"

shift 3

EOS_BASE="root://junoeos01.ihep.ac.cn:1094/"
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

input_file=$(sed -n "$((PROC_ID + 1))p" "$OUTPUT_DIR/$LIST_FILE")
if [[ -z "$input_file" ]]; then
    echo "No input file found for PROC_ID=$PROC_ID"
    exit 1
fi

input_filename=$(basename "$input_file")
local_input_file="$TEMP/$input_filename"
local_vertex_file="$TEMP/${input_filename/.rtraw/.vertex.rec}"
local_track_file="$TEMP/${input_filename/.rtraw/.track.rec}"
local_output_file="$TEMP/${input_filename/.rtraw/.output.root}"
output_file="$OUTPUT_DIR/$(basename "$local_output_file")"

echo "Output filename: $output_file"

if [[ "$SKIP_IF_EXISTS" == true && -f "$output_file" ]]; then
    echo "Skipping PROC_ID=$PROC_ID — output already exists: $output_file"
    exit 0
fi

echo "Copying file from EOS: $input_file"
xrdcp "${EOS_BASE}${input_file}" "$local_input_file"

(
    source /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.5.0/setup.sh
    echo "TUTORIALROOT = ${TUTORIALROOT}"

    echo "Running share/tut_rtraw2rec.py for file: $local_input_file"

    python ${TUTORIALROOT}/share/tut_rtraw2rec.py \
        --loglevel Info \
        --evtmax -1 \
        --method qctr \
        --global-tag MixedPhase_J25.7.2 \
        --waverec-method cotiwaverec \
        --Calib 1 \
        --pmtcalibsvc-ChargeAlgType 0 \
        --pmtcalibsvc-ReadDB 1 \
        --pmtcalibsvc-DBcur 20250210 \
        --input "$local_input_file" \
        --output "$local_vertex_file" \
        --output-stream /Event/CdLpmtCalib:on \
        --output-stream /Event/CdTrigger:on \
        --output-stream /Event/WpCalib:on \
        --output-stream /Event/WpTrigger:on
)

(
    source /afs/ihep.ac.cn/users/t/traymond/J25.3.0/git_junosw_J25_load.sh
    echo "TUTORIALROOT = ${TUTORIALROOT}"

    echo "Running run_muon.py for file: $local_vertex_file"

    python run_muon.py \
      --input "$local_vertex_file" \
      --output "$local_track_file" \
      "${EXTRA_ARGS[@]}"

    echo "Running run_analysis.py for file: $local_track_file"

    python run_analysis.py \
      --input "$local_track_file" \
      --output "$local_output_file" \
      "${EXTRA_ARGS[@]}"
)

echo "Copying result to $output_file"
cp "$local_output_file" "$output_file"

echo "Done (PROC_ID=$PROC_ID)"