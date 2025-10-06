#ifndef ANALYSISGROUPC_ANALYSISGROUPC_HPP_
#define ANALYSISGROUPC_ANALYSISGROUPC_HPP_

#include "SniperKernel/AlgBase.h"

#include <deque>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

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

    // Loader

    std::shared_ptr<Loader> m_loader;

    // Reconstruction tool

    std::string m_recToolName; 
    IRecMuonTool* m_recTool; 

    Params m_params; // set of parameters' key/value
    PmtTable m_pmtTable;

    // Properties

    double m_sigmaPmt20inch;
    double m_sigmaPmt3inch;
    double m_sigmaPmtTt;

    bool m_flagUse20inch;
    bool m_flagUse3inch;
    int m_chosenDetectors;

    bool m_useJointLoader;
    std::pair<double, double> m_loaderTimeWindow;

    // Muon selection variable

    bool m_reconstruct_muon_mode;
    bool m_first_reconstruction_file;
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

    ContextFileTracker m_contextTracker;
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