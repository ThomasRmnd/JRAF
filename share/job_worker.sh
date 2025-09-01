#!/bin/bash
set -e

source /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.5.0/setup.sh

# Arguments: <run_number> <proc_id> <start_file> [other args...]
run_number=$1
proc_id=$2
start_file=$3

source /afs/ihep.ac.cn/users/t/traymond/J25.3.0/git_junosw_J25_load.sh