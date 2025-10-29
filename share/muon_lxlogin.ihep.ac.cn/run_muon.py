import argparse
import json
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="Input filepath")
parser.add_argument("--input-rtraw", type=str, help="Input RTRAW filepath")
parser.add_argument("--output", type=str, help="Output filepath")
parser.add_argument("--property-file", type=str, help="Filepath of the property file")
parser.add_argument("--tt-reco-filepath", type=str, help="TT reco filepath")
parser.add_argument("--time-window", nargs=2, type=float, metavar=("START", "END"), help="Buffer time window")
parser.add_argument("--log-level", type=int, default=1, help="Log level (default: 1)")
args = parser.parse_args()

ifilepath = args.input
rtraw_ifilepath = args.input_rtraw
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

# ~~~~~~~~~~ PMTParamSvc ~~~~~~~~~~
pmt_param_svc = task.createSvc("PMTParamSvc")

# ~~~~~~~~~~ TTGeomSvc ~~~~~~~~~~
tt_geom_svc = task.createSvc("TTGeomSvc")

# ~~~~~~~~~~ RootIOSvc ~~~~~~~~~~
import RootIOSvc
ri_svc = task.createSvc("RootInputSvc/InputSvc")
ri_svc.property("InputFile").set(ifilepath)
ri_svc.property("InputCorrelationFile").set(rtraw_ifilepath)

output_streams = {
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
import WpMuonClassifyRecTool
alg = AnalysisGroupC.createAlg(task)
alg.setLogLevel(loglevel)

alg.useRecTool("CdWpTtChi2RecTool")
alg.useClassifyTool("WpMuonClassifyRecTool")

alg.useLoader("JointLoader")
alg.loader.property("TimeWindow").set([-500.0, 500.0]) # ns

alg.useCdFiller("CdRangeFiller")
alg.cdfiller.property("Pmt3inchTimeReso").set(15.0)
alg.cdfiller.property("Pmt20inchTimeReso").set(8.0)

alg.useWpFiller("WpRangeFiller")
alg.wpfiller.property("PmtTimeReso").set(8.0)

alg.useClassifyLoader()
alg.classify_wpfiller.property("PmtTimeReso").set(8.0)

alg.property("TtRecoFilepath").set(args.tt_reco_filepath)
alg.property("ReconstructMuonMode").set(True)

if args.property_file:
    try:
        with open(args.property_file, "r") as f:
            props = json.load(f)
        for key, value in props.items():
            print(f"[INFO] Setting alg.rectool.property('{key}') = {value}")
            alg.rectool.property(key).set(value)
    except Exception as e:
        print(f"[ERROR] Failed to load or parse property file: {e}", file=sys.stderr)
        sys.exit(1)

alg.classifytool.property("WpMuonClassifyRecToolInitialChargeCut").set(28.0)
alg.classifytool.property("WpMuonClassifyRecToolMaxChargeThreshold").set(200.0)
alg.classifytool.property("WpMuonClassifyRecToolDistanceThreshold").set(6500.0)
alg.classifytool.property("UseAdditionalGainCorrection").set(True)
alg.classifytool.property("AdditionalGainCorrectionPath").set("/junofs/users/traymond/data/WpClassifyMuonRecTool/RatioCopyNo.txt")

task.setEvtMax(-1)
if not task.run():
    print("Task ran failed!", flush=True)
    sys.exit(1)
sys.exit(0)