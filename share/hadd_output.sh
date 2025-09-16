#!/bin/bash

usage() {
  echo "Usage:"
  echo "  $0 -f PhysicsRun.txt -s 20250826_20250907"
  echo "  $0 -f PhysicsRun.txt -s 20250826_20250907 -r <min_run> -R <max_run>"
  echo
  echo "Options:"
  echo "  -f <file>      PhysicsRun.txt (run_number and output_path per line)"
  echo "  -s <summary>   Summary directory name (e.g. 20250826_20250907)"
  echo "  -r <run_min>   Minimum run number (inclusive)"
  echo "  -R <run_max>   Maximum run number (inclusive)"
  exit 1
}

# --- Parse args ---
while getopts "f:s:r:R:" opt; do
  case $opt in
    f) run_file=$OPTARG ;;
    s) summary_dir=$OPTARG ;;
    r) run_min=$OPTARG ;;
    R) run_max=$OPTARG ;;
    *) usage ;;
  esac
done

[ -z "$run_file" ] && usage
[ -z "$summary_dir" ] && usage

summary_base="/junofs/users/traymond/analysis/ibd/summary/$summary_dir"
mkdir -p "$summary_base"

# --- Loop over PhysicsRun.txt ---
while read -r run_number outpath; do
  # Skip empty/comment lines
  [[ -z "$run_number" || "$run_number" =~ ^# ]] && continue

  # If run range is specified, filter
  if [ -n "$run_min" ] && [ -n "$run_max" ]; then
    if (( run_number < run_min || run_number > run_max )); then
      continue
    fi
  fi

  echo "Processing run $run_number from $outpath..."

  # Find all matching ROOT files
  files=$(ls ${outpath}/RUN.${run_number}.JUNODAQ.Physics.ds-2.global_trigger.*.output.root 2>/dev/null)
  if [ -z "$files" ]; then
    echo "  No files found for run $run_number"
    continue
  fi

  # Extract YYYYMMDD from the first file
  first_file=$(echo $files | awk '{print $1}')
  yyyymmdd=$(basename "$first_file" | cut -d'.' -f7 | cut -c1-8)

  # Extract version string (between last '_' and '.output.root')
  version=$(basename "$first_file" | sed -E 's/.*_([^_]+_[^_]+)\.output\.root/\1/')

  # Construct output filename (drop HHMMSS and file index, keep version)
  output_file="${summary_base}/RUN.${run_number}.JUNODAQ.Physics.ds-2.global_trigger.${yyyymmdd}.${version}.output.root"

  echo "  Merging -> $output_file"
  hadd -f "$output_file" $files

done < "$run_file"
