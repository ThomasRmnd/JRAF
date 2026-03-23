#ifndef JRAF_UTILS_TRACKSAVER_HPP_
#define JRAF_UTILS_TRACKSAVER_HPP_

#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>

#include "SniperKernel/SniperLog.h"

#include "EvtNavigator/EvtNavigator.h"

#include "event/Event.hpp"
#include "event/EventCache.hpp"

struct TrackSaver {

    std::string filename;
    std::string treename = "muons";
    TFile* file = nullptr;
    TTree* tree = nullptr;

    int run_id = 0;
    time_t sec = 0l;
    int nsec = 0;
    double totq_cd = 0.0;
    double totq_wp = 0.0;

    std::vector<std::string> method;
    std::vector<unsigned char> det;
    std::vector<double> quality;

    std::vector<double> iposx;
    std::vector<double> iposy;
    std::vector<double> iposz;
    std::vector<double> fposx;
    std::vector<double> fposy;
    std::vector<double> fposz;


    bool init() {
        file = TFile::Open(filename.c_str(), "RECREATE");
        if (!file) {
            LogWarn << "Cannot open ROOT file " << filename << ". Skipping track saving\n";
            return true;
        }
        tree = new TTree(treename.c_str(), treename.c_str());
        if (!tree) {
            LogError << "Cannot create track TTree\n";
            return false;
        }
        
        tree->Branch("run_id", &run_id);
        tree->Branch("sec", &sec);
        tree->Branch("nsec", &nsec);
        tree->Branch("totq_cd", &totq_cd);
        tree->Branch("totq_wp", &totq_wp);
        
        tree->Branch("method", &method);
        tree->Branch("det", &det);
        tree->Branch("quality", &quality);
        
        tree->Branch("iposx", &iposx);
        tree->Branch("iposy", &iposy);
        tree->Branch("iposz", &iposz);
        tree->Branch("fposx", &fposx);
        tree->Branch("fposy", &fposy);
        tree->Branch("fposz", &fposz);
        
        return true;
    }

    void reset() {
        run_id = 0;
        sec = 0l;
        nsec = 0;
        totq_cd = 0.0;
        totq_wp = 0.0;
        method.clear();
        det.clear();
        quality.clear();
        iposx.clear();
        iposy.clear();
        iposz.clear();
        fposx.clear();
        fposy.clear();
        fposz.clear();
    }

    void fill(JM::EvtNavigator* nav) {
        if (!tree) return;
        std::shared_ptr<Event> evt = EventCache::load(nav);
        if (evt->tracks.empty()) return;
        reset();
        run_id = evt->run_id;
        sec = evt->ts.GetSec();
        nsec = evt->ts.GetNanoSec();
        totq_cd = evt->totq_cd;
        totq_wp = evt->totq_wp;
        for (const track& t : evt->tracks) {
            method.push_back(t.method);
            det.push_back(static_cast<unsigned char>(t.det));
            quality.push_back(t.quality);
            iposx.push_back(t.ipos.x);
            iposy.push_back(t.ipos.y);
            iposz.push_back(t.ipos.z);
            fposx.push_back(t.fpos.x);
            fposy.push_back(t.fpos.y);
            fposz.push_back(t.fpos.z);
        }
        if (method.empty()) return;
        tree->Fill();
    }

    bool save() {
        if (!file || !tree) return true;
        file->cd();
        tree->Write();
        file->Close();
        return true;
    }

};

#endif // JRAF_UTILS_TRACKSAVER_HPP_