#ifndef JRAF_VETO_VETOTIMESAVER_HPP_
#define JRAF_VETO_VETOTIMESAVER_HPP_

#include <TTree.h>

#include "SniperKernel/SniperLog.h"

#include "Context/TimeStamp.h"

#include "veto/Veto.hpp"

struct VetoWindow {
    TimeStamp start;
    TimeStamp end;
    VetoType type;
};

struct VetoTimeSaver {

    std::unordered_map<VetoType, TimeStamp> veto_map {
        {VetoType::BeginningOfJob, TimeStamp{1, 200000000}},
        {VetoType::MissingHeaders, TimeStamp{0, 5000000}},
        {VetoType::BigGaps, TimeStamp{1, 200000000}},
        {VetoType::Muon, TimeStamp{0, 5000000}}
    };

    std::deque<VetoWindow> active_vetoes;

    TTree* tree = nullptr;

    int run_id = 0;
    time_t sec = 0l;
    int nsec = 0;
    unsigned char veto_type = 0;

    bool init() {
        tree = new TTree("Veto", "Veto");
        if (!tree) {
            LogError << "Cannot create veto TTree\n";
            return false;
        }
        tree->Branch("run_id", &run_id);
        tree->Branch("sec", &sec);
        tree->Branch("nsec", &nsec);
        tree->Branch("veto_type", &veto_type);
        return true;
    }

    bool create(const TimeStamp& ts, VetoType type, int run) {
        std::unordered_map<VetoType, TimeStamp>::const_iterator it = veto_map.find(type);
        if (it == veto_map.end()) {
            return false;
        }

        TimeStamp duration = it->second;
        VetoWindow window{ts, ts + duration, type};
        active_vetoes.push_back(window);

        run_id = run;
        sec = ts.GetSec();
        nsec = ts.GetNanoSec();
        veto_type = static_cast<unsigned char>(it->first);
        tree->Fill();
        
        return true;
    }

    bool create_no_veto(const TimeStamp& ts, VetoType type, int run) {
        std::unordered_map<VetoType, TimeStamp>::const_iterator it = veto_map.find(type);
        if (it == veto_map.end()) {
            return false;
        }
        
        run_id = run;
        sec = ts.GetSec();
        nsec = ts.GetNanoSec();
        veto_type = static_cast<unsigned char>(it->first);
        tree->Fill();
        
        return true;
    }

    bool inVeto(const TimeStamp& ts) {
        cleanupExpired(ts);
        for (const VetoWindow& veto : active_vetoes) {
            if (veto.start <= ts && ts < veto.end) {
                return true;
            }
        }
        return false;
    }

    void cleanupExpired(const TimeStamp& ts) {
        while (!active_vetoes.empty() && ts >= active_vetoes.front().end) {
            active_vetoes.pop_front();
        }
    }

    bool write() {
        if (!tree) return false;
        tree->Write();
        return true;
    }

};

#endif // JRAF_VETO_VETOTIMESAVER_HPP_