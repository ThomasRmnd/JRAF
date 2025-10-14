#!/bin/bash
set -e

PROC_ID="$1"
RANGE_START="$2"
RANGE_END="$3"
RUN_NUMBER="$4"
LIST_BASE="$5"
shift 5
EXTRA_ARGS=("$@")

FILE_NUMBER=$((RANGE_START + PROC_ID))

if (( FILE_NUMBER > RANGE_END )); then
    echo "ERROR: FILE_NUMBER=$FILE_NUMBER exceeds RANGE_END=$RANGE_END"
    exit 1
fi

EOS_BASE="root://junoeos01.ihep.ac.cn/"
RTRAW_LIST_FILE="$LIST_BASE/rtraw_list/run_${RUN_NUMBER}.txt"
ESD_LIST_FILE="$LIST_BASE/esd_list/run_${RUN_NUMBER}.txt"

mapfile -t rtraw_list < <(xrdfs "$EOS_BASE" cat "$RTRAW_LIST_FILE")
mapfile -t esd_list < <(xrdfs "$EOS_BASE" cat "$ESD_LIST_FILE")

if (( ${#rtraw_list[@]} != ${#esd_list[@]} )); then
    echo "Error: mismatch in rtraw (${#rtraw_list[@]}) vs esd (${#esd_list[@]})"
    exit 1
fi

input_rtraw_file=""
input_esd_file=""

for ((i=0; i<${#esd_list[@]}; i++)); do
    esd_file="${esd_list[$i]}"
    fname=${esd_file##*/}
    
    if [[ $fname =~ \.[0-9]{14}\.([0-9]+)_ ]]; then
        file_number=${BASH_REMATCH[1]}
        file_number=$((10#$file_number))  # strip leading zeros
        
        if (( file_number == FILE_NUMBER )); then
            input_esd_file="$esd_file"
            input_rtraw_file="${rtraw_list[$i]}"
            break
        fi
    fi
done

if [[ -z "$input_rtraw_file" || -z "$input_esd_file" ]]; then
    echo "No input file found for FILE_NUMBER=$FILE_NUMBER"
    exit 1
fi

input_rtraw_filename=$(basename "$input_rtraw_file")
input_esd_filename=$(basename "$input_esd_file")

output_path=""
if [[ "$input_esd_file" =~ /eos/juno/esd/([^/]+)/([0-9]{4})/([0-9]{4})/RUN\.([0-9]+)\. ]]; then
    esd_version="${BASH_REMATCH[1]}"
    year="${BASH_REMATCH[2]}"
    monthday="${BASH_REMATCH[3]}"
    output_path="/scratchfs/juno/traymond/analysis/ibd/${year}/${monthday}"
elif [[ "$input_esd_file" =~ /eos/juno-kup/([^/]+)/([^/]+)/([^/]+)/([^/]+)/([^/]+)/RUN\.([0-9]+)\. ]]; then
    campaign="${BASH_REMATCH[1]}"
    stream="${BASH_REMATCH[2]}"
    run_bucket="${BASH_REMATCH[3]}"
    run_group="${BASH_REMATCH[4]}"
    run_number="${BASH_REMATCH[5]}"
    output_path="/scratchfs/juno/traymond/analysis/ibd/${run_bucket}/${run_group}/${RUN_NUMBER}"
else
    echo "Error: unrecognized esd file path format: $input_esd_file"
    exit 1
fi
mkdir -p "$output_path"

echo "output_path: ${output_path}" 1>&2

# input esd filename should have the format: 
# - /eos/juno/esd/<esd version>/<year>/<month><day>/RUN.<RUN_NUMBER>.*
# - /eos/juno-kup/<campaign>/<stream>/<run_bucket>/<run_group>/<run>/RUN.<RUN_NUMBER>.*

# we want for tt_reco_filepath:
# - /eos/juno/dirac/juno/user/j/jpandre_1/tt_data_auto/<year>/<month><day>/RUN.<RUN_NUMBER>.*
# - /eos/juno/dirac/juno/user/j/jpandre_1/tt_data_auto/<run_bucket>/<run_group>/<run>/RUN.<RUN_NUMBER>.*

tt_reco_filepath=""
if [[ "$input_esd_file" =~ /eos/juno/esd/([^/]+)/([0-9]{4})/([0-9]{4})/RUN\.([0-9]+)\. ]]; then
    esd_version="${BASH_REMATCH[1]}"
    year="${BASH_REMATCH[2]}"
    monthday="${BASH_REMATCH[3]}"
    tt_reco_filepath="${EOS_BASE}/eos/juno/dirac/juno/user/j/jpandre_1/tt_data_auto/${year}/${monthday}/RUN.${RUN_NUMBER}.*"
elif [[ "$input_esd_file" =~ /eos/juno-kup/([^/]+)/([^/]+)/([^/]+)/([^/]+)/([^/]+)/RUN\.([0-9]+)\. ]]; then
    campaign="${BASH_REMATCH[1]}"
    stream="${BASH_REMATCH[2]}"
    run_bucket="${BASH_REMATCH[3]}"
    run_group="${BASH_REMATCH[4]}"
    run_number="${BASH_REMATCH[5]}"
    tt_reco_filepath="${EOS_BASE}/eos/juno/dirac/juno/user/j/jpandre_1/tt_data_auto/${run_bucket}/${run_group}/${RUN_NUMBER}/RUN.${RUN_NUMBER}.*"
else
    echo "Error: unrecognized esd file path format: $input_esd_file"
    exit 1
fi

echo "tt_reco_filepath: $tt_reco_filepath" 1>&2

local_input_rtraw_file="$TEMP/$input_rtraw_filename"
local_input_esd_file="$TEMP/$input_esd_filename"
local_output_file="$TEMP/${input_esd_filename/.esd/.track.rec}"
output_file="$output_path/$(basename "$local_output_file")"

echo "Output filename: $output_file"

echo "Copying rtraw file from EOS: $local_input_rtraw_file"
xrdcp "${input_rtraw_file}" "$local_input_rtraw_file"

echo "Copying esd file from EOS: $local_input_esd_file"
xrdcp "${input_esd_file}" "$local_input_esd_file"

source /afs/ihep.ac.cn/users/t/traymond/J25.3.0/git_junosw_J25_load.sh
echo "TUTORIALROOT = ${TUTORIALROOT}"

echo "Running run_muon.py for file: $local_input_esd_file (FILE_NUMBER=$FILE_NUMBER, PROC_ID=$PROC_ID)"

python run_muon.py \
    --input "$local_input_esd_file" \
    --input-rtraw "$local_input_rtraw_file" \
    --output "$local_output_file" \
    --tt-reco-filepath "$tt_reco_filepath" \
    "${EXTRA_ARGS[@]}"

echo "Copying result to $output_file"
cp "$local_output_file" "$output_file"

echo "Done (FILE_NUMBER=$FILE_NUMBER, PROC_ID=$PROC_ID)"