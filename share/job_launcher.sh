#!/bin/bash

EOS_BASE="root://junoeos01.ihep.ac.cn:1094/"
LIST_FILE="edm_file_list.txt"

log_level=3
time_window=("-2.0" "2.0")
skip_if_exists=false

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --input-path)
            input_path="$2"
            shift 2
            ;;
        --run-number)
            run_number="$2"
            shift 2
            ;;
        --output-path)
            output_path="$2"
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

if [[ -z "$input_path" || -z "$run_number" || -z "$output_path" ]]; then
    echo "Usage: $0 --input-path <path> --run-number <number> --output-path <path> [--file-offset <num>] [--file-range <num>] [--time-window <num> <num>]"
    exit 1
fi

echo "Listing ROOT files from EOS..."
mapfile -t file_list < <(xrdfs "$EOS_BASE" ls "$input_path")

echo "Number of file before applying run number: ${#file_list[@]}"

file_list=($(printf "%s\n" "${file_list[@]}" | grep "RUN\.${run_number}\.*\.esd"))

echo "Number of file after applying run number: ${#file_list[@]}"

if [[ -z "$file_offset" ]]; then
    file_offset=0
fi

if [[ -z "$file_range" ]]; then
    file_range=$(( ${#file_list[@]} - file_offset ))
fi

if ! [[ "$file_offset" =~ ^[0-9]+$ && "$file_range" =~ ^[0-9]+$ ]]; then
    echo "file-offset and file-range must be non-negative integers"
    exit 1
fi

echo "Total number of file: [$file_offset, $file_range]"

file_list=("${file_list[@]:$file_offset:$file_range}")
job_count=${#file_list[@]}

if (( job_count == 0 )); then
    echo "No ROOT files found in $input_path"
    exit 1
fi

printf "%s\n" "${file_list[@]}" > "$output_path/$LIST_FILE"

extra_args=""

extra_args+=" --time-window ${time_window[0]} ${time_window[1]}"
extra_args+=" --log-level $log_level"

if [[ "$skip_if_exists" == true ]]; then
    extra_args+=" --skip-if-exists"
fi


# if [[ -z "$property_file" ]]; then
#     property_file="/junofs/users/traymond/reconstruction/esd/properties/RUN.${run_number}.Properties.json"
# fi

# extra_args="--property-file $property_file"


# --- Submit batch jobs ---
echo "Submitting $job_count jobs with hep_sub..."
hep_sub job_worker.sh \
  -argu "%{ProcId} $LIST_FILE $output_path $extra_args" \
  -n "$job_count" \
  -cpu 1 \
  -m 4096 \
  -o "/scratchfs/juno/traymond/agrpc_${run_number}_%{ProcId}.log" \
  -e "/scratchfs/juno/traymond/agrpc_${run_number}_%{ProcId}.err" \
  -wt short \
  -name agrpc_${run_number}_batch