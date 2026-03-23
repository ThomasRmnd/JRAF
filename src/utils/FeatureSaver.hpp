#ifndef JRAF_UTILS_FEATURESAVER_HPP_
#define JRAF_UTILS_FEATURESAVER_HPP_

#include <string>

#include <TFile.h>
#include <TTree.h>

#include "SniperKernel/SniperLog.h"

struct FeatureSaver {

    std::string filename;
    std::string treename = "Features";
    TFile* file = nullptr;
    TTree* tree = nullptr;

    int run_id;
    time_t sec;
    int nsec;

    std::vector<double> iposx;
    std::vector<double> iposy;
    std::vector<double> iposz;
    std::vector<double> fposx;
    std::vector<double> fposy;
    std::vector<double> fposz;
    std::vector<double> chi2;
    std::vector<unsigned char> det; // 1 = CdWpTtChi2, 2 = CdClassify, 4 = TT

    std::vector<unsigned int> id;
    std::vector<double> fht;
    std::vector<double> totq;
    std::vector<double> q;
    std::vector<int> nhit;

    std::vector<double> pointx;
    std::vector<double> pointy;
    std::vector<double> pointz;

    bool init() {
        file = TFile::Open(filename.c_str(), "RECREATE");
        if (!file) {
            LogWarn << "Cannot open ROOT file " << filename << ". Skipping feature saving\n";
            return true;
        }
        tree = new TTree(treename.c_str(), treename.c_str());

        tree->Branch("run_id", &run_id);
        tree->Branch("sec", &sec);
        tree->Branch("nsec", &nsec);
        
        tree->Branch("iposx", &iposx);
        tree->Branch("iposy", &iposy);
        tree->Branch("iposz", &iposz);
        tree->Branch("fposx", &fposx);
        tree->Branch("fposy", &fposy);
        tree->Branch("fposz", &fposz);
        tree->Branch("chi2", &chi2);
        tree->Branch("det", &det);

        tree->Branch("id", &id);
        tree->Branch("fht", &fht);
        tree->Branch("totq", &totq);
        tree->Branch("q", &q);
        tree->Branch("nhit", &nhit);

        tree->Branch("pointx", &pointx);
        tree->Branch("pointy", &pointy);
        tree->Branch("pointz", &pointz);

        return true;  
    }
    
    void reset() {
        run_id = 0;
        sec = 0l;
        nsec = 0;

        iposx.clear();
        iposy.clear();
        iposz.clear();
        fposx.clear();
        fposy.clear();
        fposz.clear();
        chi2.clear();
        det.clear();

        id.clear();
        fht.clear();
        totq.clear();
        q.clear();
        nhit.clear();

        pointx.clear();
        pointy.clear();
        pointz.clear();
    }

    void fill() {
        if (tree && !chi2.empty() && !id.empty()) tree->Fill();
    }

    bool save() {
        if (!file || !tree) return true;
        file->cd();
        tree->Write();
        file->Close();
        return true;
    }

};

#endif // JRAF_UTILS_FEATURESAVER_HPP_