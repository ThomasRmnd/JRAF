import argparse
import json
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("path", type=str, help="Input file path")
args = parser.parse_args()

ipath = args.path
opath = args.path

if not os.path.exists(ipath):
    print(f"[ERROR] Input file does not exist: {ipath}", file=sys.stderr, flush=True)
    sys.exit(1)

if os.path.exists(opath):
    print(f"[WARN] Output file already exist: {opath}", flush=True)

print(f"[INFO] Processing file: {ipath}", flush=True)
print(f"[INFO] Output file will be: {opath}")

# === Sniper ====
import Sniper
Sniper.setLogLevel(1)
task = Sniper.TopTask("task")
task.setLogLevel(1)

# === Profiling ===
import SniperProfiling
prof = task.createSvc("SniperProfiling")
prof.setLogLevel(1)

# === BufferMemMgr ===
import BufferMemMgr
buf_mgr = task.createSvc("BufferMemMgr")
buf_mgr.property("TimeWindow").set([-1.0e-6, 1.0e-6])

# === RootIOSvc ===
import RootIOSvc
input_files = [ipath]
ri_svc = task.createSvc("RootInputSvc/InputSvc")
ri_svc.property("InputFile").set(input_files)
output_files = {
    # === Calib ===
    "/Event/CdLpmtCalib": opath,
    "/Event/CdSpmtCalib": opath,
    "/Event/WpCalib": opath,
    "/Event/TtCalib": opath,
}
ro_svc = task.createSvc("RootOutputSvc/OutputSvc")
ro_svc.property("OutputStreams").set(output_files)

task.setEvtMax(-1)
# task.show()
if (task.run()):
    print("Task finished successfully!")
    sys.exit(0)
else:
    print("Task failed!")
    sys.exit(1)