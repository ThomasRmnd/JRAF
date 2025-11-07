import numpy as np
import uproot

def sum_daq_time(filename):
    """Special analysis for Thomas' DAQTree"""
    with uproot.open(filename) as f:
        if "DAQTree" not in f:
            print("[INFO] No DAQTree found in file")
            return
        tree = f["DAQTree"]

        daq_sec = tree["daq_sec"].array(library="np")
        daq_nsec = tree["daq_nsec"].array(library="np")
        daq_ts = daq_sec.astype(np.float64) + daq_nsec * 1e-9
        daq_sum = daq_ts.sum()

        daq_days = daq_sum / 86400

        return daq_days

def sum_daq_muon(filename):
    """Special analysis for Thomas' DAQTree"""
    with uproot.open(filename) as f:
        if "DAQTree" not in f:
            print("[INFO] No DAQTree found in file")
            return
        tree = f["DAQTree"]

        muveto_sec = tree["muveto_sec"].array(library="np")
        muveto_nsec = tree["muveto_nsec"].array(library="np")
        muveto_ts = muveto_sec.astype(np.float64) + muveto_nsec * 1e-9
        muveto_sum = muveto_ts.sum()

        muveto_days = muveto_sum / 86400

        return muveto_days
