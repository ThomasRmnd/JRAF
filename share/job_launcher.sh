#!/bin/bash

EOS_BASE="root://junoeos01.ihep.ac.cn/"

log_level=3
time_window=("-2.0" "2.0")
skip_if_exists=false
list_base="/eos/juno/groups/DataQuality/P25A/Physics/goodrunlist_v2.1"

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --run-number)
            run_number="$2"
            shift 2
            ;;
        --output-path)
            output_path="$2"
            shift 2
            ;;
        --list-base)
            list_base="$2"
            shift 2
            ;;
        --file-offset)
            file_offset="$2"
            shift 2
            ;;
        --file-range)
            file_range="$2"
            shift 2
            ;;
        --time-window)
            time_window=("$2" "$3")
            shift 3
            ;;
        --log-level)
            log_level="$2"
            shift 2
            ;;
        --property-file)
            property_file="$2"
            shift 2
            ;;
        --skip-if-exists)
            skip_if_exists=true
            shift
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

if [[ -z "$run_number" || -z "$output_path" ]]; then
    echo "Usage: $0 --run-number <number> --output-path <path> [--file-offset <num>] [--file-range <num>] [--time-window <num> <num>]"
    exit 1
fi

RTRAW_LIST_FILE="${list_base}/rtraw_list/run_${run_number}.txt"
ESD_LIST_FILE="${list_base}/esd_list/run_${run_number}.txt"

echo "Listing ROOT files from EOS..."
mapfile -t rtraw_list < <(xrdfs "$EOS_BASE" cat "$RTRAW_LIST_FILE")
mapfile -t esd_list   < <(xrdfs "$EOS_BASE" cat "$ESD_LIST_FILE")

echo "Number of rtraw file: ${#rtraw_list[@]}"
echo "Number of esd file: ${#esd_list[@]}"

if [[ -z "$file_offset" ]]; then
    file_offset=0
fi

if [[ -z "$file_range" ]]; then
    file_range=$(( ${#rtraw_list[@]} - file_offset ))
fi

if ! [[ "$file_offset" =~ ^[0-9]+$ && "$file_range" =~ ^[0-9]+$ ]]; then
    echo "file-offset and file-range must be non-negative integers"
    exit 1
fi

echo "Total number of file: [$file_offset, $file_range]"

rtraw_list=("${rtraw_list[@]:$file_offset:$file_range}")
esd_list=("${esd_list[@]:$file_offset:$file_range}")

job_count_rtraw=${#rtraw_list[@]}
job_count_esd=${#esd_list[@]}

if (( job_count_rtraw != job_count_esd )); then
    echo "Error: mismatch between rtraw files ($job_count_rtraw) and esd files ($job_count_esd)"
    exit 1
fi

job_count=${#rtraw_list[@]}

if (( job_count == 0 )); then
    echo "No ROOT files found in $input_path"
    exit 1
fi

extra_args=" --time-window ${time_window[0]} ${time_window[1]} --log-level $log_level"

if [[ "$skip_if_exists" == true ]]; then
    extra_args+=" --skip-if-exists"
fi


# if [[ -z "$property_file" ]]; then
#     property_file="/junofs/users/traymond/reconstruction/esd/properties/RUN.${run_number}.Properties.json"
# fi

# extra_args="--property-file $property_file"

mkdir -p "$output_path"

# --- Submit batch jobs ---
echo "Submitting $job_count jobs with hep_sub..."
hep_sub job_worker.sh \
  -argu "%{ProcId} $run_number $list_base $output_path $extra_args" \
  -n "$job_count" \
  -cpu 1 \
  -m 4096 \
  -o "/scratchfs/juno/traymond/agrpc_${run_number}_%{ProcId}.log" \
  -e "/scratchfs/juno/traymond/agrpc_${run_number}_%{ProcId}.err" \
  -name agrpc_${run_number}_batch
#   -wt short \
#   -o "/scratchfs/juno/traymond/agrpc_${run_number}_%{ProcId}.log" \
#   -e "/scratchfs/juno/traymond/agrpc_${run_number}_%{ProcId}.err" \