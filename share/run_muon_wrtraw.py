import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="Input filepath")
parser.add_argument("--input-rtraw", type=str, help="Input RTRAW filepath")
parser.add_argument("--output", type=str, help="Output filepath")
parser.add_argument("--time-window", nargs=2, type=float, metavar=("START", "END"), help="Buffer time window")
parser.add_argument("--log-level", type=int, default=1, help="Log level (default: 1)")
args = parser.parse_args()

ifilepath = args.input
ofilepath = args.output
lower_tw, upper_tw = args.time_window
loglevel = args.log_level

# ~~~~~~~~~~ Sniper ~~~~~~~~~~
import Sniper
Sniper.setLogLevel(loglevel)
task = Sniper.TopTask("task")
task.setLogLevel(loglevel)

# ~~~~~~~~~~ Profiling ~~~~~~~~~~
import SniperProfiling
prof = task.createSvc("SniperProfiling")
prof.setLogLevel(loglevel)

# ~~~~~~~~~~ BufferMemMgr ~~~~~~~~~~
import BufferMemMgr
buf_mgr = task.createSvc("BufferMemMgr")
buf_mgr.property("TimeWindow").set([lower_tw, upper_tw])

# ~~~~~~~~~~ Geometry ~~~~~~~~~~ 
import Geometry
geom = task.createSvc("RecGeomSvc")
geom.property("GeomFile").set("default")
geom.property("GeomPathInRoot").set("JunoGeom")
geom.property("FastInit").set(True)

# ~~~~~~~~~~ TTGeomSvc ~~~~~~~~~~
tt_geom_svc = task.createSvc("TTGeomSvc")

# ~~~~~~~~~~ RootIOSvc ~~~~~~~~~~
import RootIOSvc
input_file = [ifilepath]
ri_svc = task.createSvc("RootInputSvc/InputSvc")
ri_svc.property("InputFile").set(input_file)
ri_svc.property("InputCorrelationFile").set(args.input_rtraw)

output_streams = {
    # === Calib ===
    "/Event/CdLpmtCalib": ofilepath,
    "/Event/CdSpmtCalib": ofilepath,
    "/Event/WpCalib": ofilepath,
    "/Event/TtCalib": ofilepath,
    # === Rec ===
    "/Event/CdVertexRec": ofilepath,
    "/Event/CdTrackRec": ofilepath,
    "/Event/WpRec": ofilepath,
    "/Event/TtRec": ofilepath
}
ro_svc = task.createSvc("RootOutputSvc/OutputSvc")
ro_svc.property("OutputStreams").set(output_streams)

# ~~~~~~~~~~ AnalysisGroupC ~~~~~~~~~~
import AnalysisGroupC
import CdWpTtChi2RecTool
alg = AnalysisGroupC.createAlg(task)
alg.setLogLevel(1)
alg.useRecTool("CdWpTtChi2RecTool")

alg.property("Pmt3inchTimeReso").set(15.0) # 15.0
alg.property("Pmt20inchTimeReso").set(8.0)
alg.property("PmtTTTimeReso").set(2.0)
alg.property("Use3inchPMT").set(True)
alg.property("Use20inchPMT").set(True)
alg.property("ChosenDetectors").set(3) # 1: CD, 2: WP, 4: TT
alg.property("UseJointLoader").set(True) 
alg.property("LoaderTimeWindow").set([-500.0, 500.0]) # ns
alg.property("ReconstructMuonMode").set(True)

task.setEvtMax(-1)
if not task.run():
    print("Task ran failed!", flush=True)
    sys.exit(1)
sys.exit(0)