import argparse
import json
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, nargs="+", help="Input filepath")
parser.add_argument("--input-correlation", type=str, nargs="+", default=[], help="Input correlation filepath")
parser.add_argument("--output", type=str, help="Output filepath")
parser.add_argument("--context-previous-filename", type=str, default="", help="Context previous filename")
parser.add_argument("--context-next-filename", type=str, default="", help="Context next filename")
parser.add_argument("--tt-reco-filepath", type=str, default="", help="TT reco filepath")
parser.add_argument("--reco-output", type=str, default="", help="Muon reconstruction output filepath")
parser.add_argument("--feature-output", type=str, default="", help="Feature output filepath")
parser.add_argument("--property-file", type=str, default="", help="Filepath of the property file")
parser.add_argument("--time-window", nargs=2, type=float, default=(-2.0, 2.0), help="Buffer time window (default=[-2.0, 2.0])")
parser.add_argument("--log-level", type=int, default=3, help="Log level (default: 3)")
parser.add_argument("--cluster", type=str, choices=["CC-IN2P3", "IHEP"], help="Cluster name currently running on")
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
pmt_param_svc = task.createSvc("PMTParamSvc")
tt_geom_svc = task.createSvc("TTGeomSvc")

# ~~~~~~~~~~ RootIOSvc ~~~~~~~~~~
import RootIOSvc
ri_svc = task.createSvc("RootInputSvc/InputSvc")
ri_svc.property("InputFile").set(ifilepath)
if args.input_correlation:
    ri_svc.property("InputCorrelationFile").set(args.input_correlation)

# ~~~~~~~~~~ JRAF ~~~~~~~~~~
import JRAF
alg = JRAF.createAlg(task)
alg.setLogLevel(loglevel)

alg.property("ContextPreviousFilename").set(args.context_previous_filename)
alg.property("ContextNextFilename").set(args.context_next_filename)
alg.property("OutputFilename").set(ofilepath)
alg.property("RecoTrackOutputFilename").set(args.reco_output)
alg.property("FeatureOutputFilename").set(args.feature_output)

alg.useEventBuilder("EventBuilder")
alg.eventbuilder.property("TtReco").set(args.tt_reco_filepath)

alg.useMuonTagger("MuonTagger")

alg.useLoader("JointLoader")
alg.loader.property("TimeWindow").set([-500.0, 500.0]) # ns

# alg.useCdFiller("CdRangeFiller")
# alg.cdfiller.property("Pmt3inchTimeReso").set(15.0)
# alg.cdfiller.property("Pmt20inchTimeReso").set(8.0)
alg.useCdFiller("CdLRangeFiller")
alg.cdfiller.property("PmtTimeReso").set(8.0)

alg.useWpFiller("WpRangeFiller")
alg.wpfiller.property("PmtTimeReso").set(8.0)

import CdWpTtChi2RecTool

alg.useRecTool("CdWpTtChi2RecTool")

alg.rectool.property("UseAutomaticFactory").set(False)
alg.rectool.property("WaterPhase").set(False)
alg.rectool.property("ManualReconstructionMode").set(1)
if args.property_file:
    alg.rectool.property("ConfigFile").set(args.property_file)


# if args.property_file:
#     try:
#         with open(args.property_file, "r") as f:
#             props = json.load(f)
#         for key, value in props.items():
#             print(f"[INFO] Setting alg.rectool.property('{key}') = {value}")
#             alg.rectool.property(key).set(value)
#     except Exception as e:
#         print(f"[ERROR] Failed to load or parse property file: {e}", file=sys.stderr)
#         sys.exit(1)

task.setEvtMax(-1)
if not task.run():
    print("Task ran failed!", flush=True)
    sys.exit(1)
sys.exit(0)