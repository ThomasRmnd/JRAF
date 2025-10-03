#!/bin/bash

EOS_BASE="root://junoeos01.ihep.ac.cn/"

list_base="/eos/juno/groups/DataQuality/P25A/Physics/goodrunlist_v2.1"
file_range=100
time_window=("-2.0" "2.0")
log_level=3

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

output_path="${output_path/\/junofs\/users/\/scratchfs\/juno}"

RTRAW_LIST_FILE="${list_base}/rtraw_list/run_${run_number}.txt"
ESD_LIST_FILE="${list_base}/esd_list/run_${run_number}.txt"

echo "Listing ROOT files from EOS..."
mapfile -t rtraw_list < <(xrdfs "$EOS_BASE" cat "$RTRAW_LIST_FILE")
mapfile -t esd_list   < <(xrdfs "$EOS_BASE" cat "$ESD_LIST_FILE")

echo "Number of rtraw file: ${#rtraw_list[@]}"
echo "Number of esd file: ${#esd_list[@]}"

if ! [[ "$file_range" =~ ^[0-9]+$ ]]; then
    echo "file-range must be a non-negative integer"
    exit 1
fi

ranges=()
range_start=""
prev_num=""
count=0

for f in "${esd_list[@]}"; do
    fname=${f##*/}
    echo "Current file: $fname"

    if [[ $fname =~ \.[0-9]{14}\.([0-9]+)_ ]]; then
        file_number=${BASH_REMATCH[1]}
        file_number=$((10#$file_number))  # strip leading zeros
    else
        echo "Warning: could not extract file number from $fname" >&2
        continue
    fi

    if [[ -z "$range_start" ]]; then
        range_start=$file_number
        count=1
    else
        if (( file_number == prev_num + 1 && count < $file_range )); then
            ((count++))
        else
            ranges+=("$range_start-$prev_num")
            range_start=$file_number
            count=1
        fi
    fi
    prev_num=$file_number
done
ranges+=("$range_start-$prev_num")

extra_args=" --time-window ${time_window[0]} ${time_window[1]} --log-level $log_level"

# if [[ -z "$property_file" ]]; then
#     property_file="/junofs/users/traymond/reconstruction/esd/properties/RUN.${run_number}.Properties.json"
# fi

# extra_args="--property-file $property_file"

mkdir -p "$output_path"

for r in "${ranges[@]}"; do
    start=${r%-*}
    end=${r#*-}
    n_jobs=$((end - start + 1))  # +1 because range is inclusive
    
    echo "Submitting $n_jobs parallel jobs for run $run_number range $start-$end"
    
    hep_sub job_worker_multifiles.sh \
        -argu "%{ProcId} $start $end $run_number $list_base $output_path $extra_args" \
        -n "$n_jobs" \
        -cpu 1 \
        -m 4096 \
        -o "/scratchfs/juno/traymond/agrpc_${run_number}_${start}_${end}_%{ProcId}.log" \
        -e "/scratchfs/juno/traymond/agrpc_${run_number}_${start}_${end}_%{ProcId}.err" \
        -name "agrpc_${run_number}_${start}_${end}_batch"
done