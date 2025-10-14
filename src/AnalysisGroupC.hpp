#ifndef ANALYSISGROUPC_ANALYSISGROUPC_HPP_
#define ANALYSISGROUPC_ANALYSISGROUPC_HPP_

#include "SniperKernel/AlgBase.h"

#include <deque>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TTree.h>

#include "Context/TimeStamp.h"
#include "EvtNavigator/NavBuffer.h"
#include "Geometry/IPMTParamSvc.h"
#include "Geometry/IRecGeomSvc.hh"
#include "Geometry/ITTGeomSvc.hh"
#include "RootIOSvc/RootInputSvc.h"
#include "RecTools/IRecMuonTool.h"
#include "RecTools/PmtProp.h"

#include "analysis/Analysis.hpp"
#include "loader/Loader.hpp"

struct ContextFileTracker {

    std::string target;
    std::string current;
    std::string next;
    bool change = false;

    bool isTarget(RootInputSvc* iptSvc) {
        if (target.empty()) return true;
        if (current.empty()) {
            current = std::filesystem::path(iptSvc->getInputStream("EvtNavigator")->streamname()).filename().string();
            next = current;
        }
        if (change) {
            current = next;
            change = false;
        }
        std::string filename = std::filesystem::path(iptSvc->getInputStream("EvtNavigator")->streamname()).filename().string();
        if (filename != current) {
            next = filename;
            change = true;
        }
        LogInfo << "Current filename: " << current << '\n';
        LogInfo << "Target filename: " << target << '\n';
        return current == target;
    }

};

struct TtRecoFile {

    std::string filename;
    std::string treename = "TT";
    TChain* chain = nullptr;
    int cur_idx = 0;
    Int_t ntracks;
    TTimeStamp* ts = nullptr;
    Int_t npts[100];
    Double_t coeff0[100], coeff1[100], coeff2[100], coeff3[100], coeff4[100], coeff5[100];
    Double_t chi2[100];

    bool init() {
        chain = new TChain(treename.c_str());
        if (!chain) {
            LogError << "Cannot create chain with name " << treename << '\n';
            return false;
        }
        chain->Add(filename.c_str());

        chain->SetBranchAddress("NTracks", &ntracks);
        chain->SetBranchAddress("NPoints", &npts);
        chain->SetBranchAddress("start_TS", &ts);
        chain->SetBranchAddress("Coeff0", &coeff0);
        chain->SetBranchAddress("Coeff1", &coeff1);
        chain->SetBranchAddress("Coeff2", &coeff2);
        chain->SetBranchAddress("Coeff3", &coeff3);
        chain->SetBranchAddress("Coeff4", &coeff4);
        chain->SetBranchAddress("Coeff5", &coeff5);
        chain->SetBranchAddress("Chi2", &chi2);

        LogInfo << "TtRecoFile has " << chain->GetEntries() << " entries\n";

        return true;
    }

    bool find(const TimeStamp& ts_) {
        long nentries = chain->GetEntries();
        TimeStamp lower_bound = ts_ - TimeStamp{0, 1000};
        TimeStamp upper_bound = ts_;
        upper_bound.Add(TimeStamp{0, 1000});
        for (; cur_idx < nentries; ++cur_idx) {
            chain->GetEntry(cur_idx);
            
            TimeStamp cur_ts{ts->GetTimeSpec()};
            if (cur_ts < lower_bound) continue;
            else if (upper_bound < cur_ts) break;

            return true;
        }
        return false;
    }

};

class AnalysisGroupC : public AlgBase {

public:

    AnalysisGroupC(const std::string& name);

    ~AnalysisGroupC() override = default;

    bool initialize() override;
    bool execute() override;
    bool finalize() override;

private:

    std::size_t m_iEvt;
    
    JM::NavBuffer* m_buf;

    // Geometry

    IRecGeomSvc* m_rgSvc;
    ITTGeomSvc* m_ttgSvc;
    IPMTParamSvc* m_pmtSvc;
    RootInputSvc* m_iptSvc;

    // Properties

    double m_sigmaPmt3inch;
    double m_sigmaPmt20inch;
    double m_sigmaPmtTt;
    bool m_flagUse3inch;
    bool m_flagUse20inch;
    int m_chosenDetectors;

    // Reconstruction

    std::string m_recToolName;
	IRecMuonTool* m_recTool;
    std::unique_ptr<Loader> m_loader;
    bool m_useJointLoader;
    std::pair<double, double> m_loaderTimeWindow;
    PmtTable m_pmtTable;
    Params m_params;

    std::string m_classifyToolName;
    IRecMuonTool* m_classifyTool;
    std::unique_ptr<Loader> m_classifyLoader;
    PmtTable m_classifyPmtTable;
    Params m_classifyParams;

    // Muon selection variable

    bool m_reconstruct_muon_mode;
    double m_cd_muon_totq_thold = 30000.0;
    double m_wp_muon_totq_thold = 400.0;
    TimeStamp m_cd_afterpulse_thold{0, 50000};
    TimeStamp m_wp_afterpulse_thold{0, 4000};
    TimeStamp m_cd_last_muon;
    TimeStamp m_wp_last_muon;

    // IBD selection variable

    TimeStamp m_muon_veto_duration = TimeStamp{0, 2000000};
    double m_prompt_lower_thold = 1500.0;
    double m_prompt_upper_thold = 20000.0;
    double m_delayed_lower_thold = 3700.0;
    double m_delayed_upper_thold = 6000.0;
    double m_distance_correlation_thold = 1500.0;
    TimeStamp m_time_correlation_lower_thold = TimeStamp{0, 1000};
    TimeStamp m_time_correlation_upper_thold = TimeStamp{0, 2000000};
    TimeStamp m_bef_multiplicity_duration = TimeStamp{0, -2000000};
    TimeStamp m_aft_multiplicity_duration = TimeStamp{0, 2000000};
    double m_fiducial_radius = 17200.0;
    double m_upper_height = 11000.0;
    double m_xyradius_thold = 3000.0;
    
    TtRecoFile m_ttRecoFile;
    ContextFileTracker m_contextTracker;
    bool m_targetIsFirst;
    TimeStamp m_targetFirstTs;
    std::vector<std::shared_ptr<Analysis>> m_analyses;

    // Output file

    std::string m_ofilename;
    TFile* m_file;

    TTree* m_daq_tree;
    time_t m_daq_sec;
    int m_daq_nsec;
    time_t m_muveto_sec;
    int m_muveto_nsec;

	bool initBufSvc(); 
	bool initGeomSvc();
    bool initRecTool();
    bool initLoader();
    bool initAnalyses();

};

#endif // ANALYSISGROUPC_ANALYSISGROUPC_HPP_