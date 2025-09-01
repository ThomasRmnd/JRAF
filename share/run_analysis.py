def run_analysis(ifilepath, ofilepath, lower_tw=2.0, upper_tw=2.0, loglevel=3):
    import os

    # ~~~~~~~~~~ Sniper ~~~~~~~~~~
    import Sniper
    Sniper.setLogLevel(loglevel)
    task = Sniper.TopTask("task")
    task.setLogLevel(loglevel)

    # ~~~~~~~~~~ Profiling ~~~~~~~~~~
    import SniperProfiling
    prof = task.createSvc("SniperProfiling/")
    prof.setLogLevel(loglevel)

    # ~~~~~~~~~~ BufferMemMgr ~~~~~~~~~~
    import BufferMemMgr
    buf_mgr = task.createSvc("BufferMemMgr/")
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

    # ~~~~~~~~~~ AnalysisGroupC ~~~~~~~~~~
    import AnalysisGroupC
    alg = AnalysisGroupC.createAlg(task)
    alg.setLogLevel(1)

    alg.property("Pmt3inchTimeReso").set(15.0) # 15.0
    alg.property("Pmt20inchTimeReso").set(8.0)
    alg.property("PmtTTTimeReso").set(2.0)
    alg.property("Use3inchPMT").set(True)
    alg.property("Use20inchPMT").set(True)
    alg.property("ChosenDetectors").set(3) # 1: CD, 2: WP, 4: TT
    alg.property("UseJointLoader").set(True) 
    alg.property("LoaderTimeWindow").set([-500.0, 500.0]) # ns
    alg.property("ReconstructionMuonMode").set(False)
    alg.property("OutputFile").set(ofilepath)

    task.setEvtMax(-1)
    if not task.run():
        print("Task ran failed!", flush=True)
        return False
    return True