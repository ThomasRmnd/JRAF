#ifndef JRAF_UTILS_DAQTIMESAVER_HPP_
#define JRAF_UTILS_DAQTIMESAVER_HPP_

#include <TTree.h>

#include "SniperKernel/SniperLog.h"

#include "Context/TimeStamp.h"

struct DAQTimeSaver {

    TTree* tree = nullptr;
    int run_id = 0;
    time_t start_sec = 0l;
    int start_nsec = 0;
    time_t duration_sec = 0l;
    int duration_nsec = 0;

    bool is_initialized = false;
    TimeStamp last_ts{0, 0};

    bool init() {
        tree = new TTree("DAQ", "DAQ");
        if (!tree) {
            LogError << "Cannot create DAQ TTree\n";
            return false;
        }
        tree->Branch("run_id", &run_id);
        tree->Branch("start_sec", &start_sec);
        tree->Branch("start_nsec", &start_nsec);
        tree->Branch("duration_sec", &duration_sec);
        tree->Branch("duration_nsec", &duration_nsec);
        return true;
    }

    bool add(const TimeStamp& ts, int run) {
        run_id = run;
        if (!is_initialized) {
            last_ts = ts;
            start_sec = ts.GetSec();
            start_nsec = ts.GetNanoSec();
            is_initialized = true;
            return true;
        }
        TimeStamp diff = ts - last_ts;
        TimeStamp daqtime{duration_sec, duration_nsec};
        daqtime.Add(diff);
        duration_sec = daqtime.GetSec();
        duration_nsec = daqtime.GetNanoSec();
        last_ts = ts;
        return true;
    };

    bool write() {
        if (!tree) return false;
        tree->Fill();
        tree->Write();
        return true;
    }

};

#endif // JRAF_UTILS_DAQTIMESAVER_HPP_