#!/bin/sh

usage() {
    cat <<EOF
Usage: $(basename "$0") --run <int> [options]

Required:
  --run <int>                  Run ID
EOF
}

if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run)           RUN="$2"; shift 2 ;;
        --help|-h)       usage; exit 0 ;;
        *) log ERROR "Unknown argument: $1"; usage; exit 1 ;;
    esac
done

if [[ -z "${RUN:-}" ]]; then
    log ERROR "--site is required {EOS|CNAF}"
    usage
    exit 1
fi

source /pbs/home/t/traymond/J25.6.1_Modified/git_junosw_J25_load.sh
python muon_rate.py "${RUN}"