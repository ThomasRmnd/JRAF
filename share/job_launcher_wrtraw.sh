#!/bin/bash

EOS_BASE="root://junoeos01.ihep.ac.cn/"

log_level=3
time_window=("-2.0" "2.0")
skip_if_exists=false

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --rtraw-path)
            rtraw_path="$2"
            shift 2
            ;;
        --esd-path)
            esd_path="$2"
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

if [[ -z "$rtraw_path" || -z "$esd_path" || -z "$run_number" || -z "$output_path" ]]; then
    echo "Usage: $0 --rtraw-path <path> --esd-path <path> --run-number <number> --output-path <path> [--file-offset <num>] [--file-range <num>] [--time-window <num> <num>]"
    exit 1
fi

RTRAW_LIST_FILE="edm_rtraw_list_${run_number}.txt"
ESD_LIST_FILE="edm_esd_list_${run_number}.txt"

echo "Listing ROOT files from EOS..."
mapfile -t rtraw_list < <(xrdfs "$EOS_BASE" ls "$rtraw_path")
mapfile -t esd_list < <(xrdfs "$EOS_BASE" ls "$esd_path")

echo "Number of rtraw file before applying run number: ${#rtraw_list[@]}"
echo "Number of esd file before applying run number: ${#esd_list[@]}"

rtraw_list=($(printf "%s\n" "${rtraw_list[@]}" | grep "RUN\.${run_number}.*\.rtraw"))
esd_list=($(printf "%s\n" "${esd_list[@]}" | grep "RUN\.${run_number}.*\.esd"))

echo "Number of rtraw file after applying run number: ${#rtraw_list[@]}"
echo "Number of esd file after applying run number: ${#esd_list[@]}"

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

printf "%s\n" "${rtraw_list[@]}" > "$output_path/$RTRAW_LIST_FILE"
printf "%s\n" "${esd_list[@]}" > "$output_path/$ESD_LIST_FILE"

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

mkdir -p "$output_path"

# --- Submit batch jobs ---
echo "Submitting $job_count jobs with hep_sub..."
hep_sub job_worker_wrtraw.sh \
  -argu "%{ProcId} $RTRAW_LIST_FILE $ESD_LIST_FILE $output_path $extra_args" \
  -n "$job_count" \
  -cpu 1 \
  -m 4096 \
  -wt short \
  -o "/scratchfs/juno/traymond/agrpc_${run_number}_%{ProcId}.log" \
  -e "/scratchfs/juno/traymond/agrpc_${run_number}_%{ProcId}.err" \
  -name agrpc_${run_number}_batch