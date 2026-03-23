#ifndef JRAF_JRAF_HPP_
#define JRAF_JRAF_HPP_

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
#include "Event/CalibPmtChannel.h"
#include "Event/CdTrackRecHeader.h"
#include "Event/CdVertexRecHeader.h"
#include "Event/OecHeader.h"
#include "Event/TtRecHeader.h"
#include "Event/WpRecHeader.h"
#include "EvtNavigator/NavBuffer.h"
#include "RecTools/IRecMuonTool.h"
#include "UtilsThomas/loader/Loader.hpp"
#include "UtilsThomas/utils/DetectorType.hpp"

#include "analysis/Analysis.hpp"
#include "event/Event.hpp"
#include "utils/ContextFileTracker.hpp"
#include "utils/DAQTimeSaver.hpp"
#include "utils/FeatureSaver.hpp"
#include "utils/TrackSaver.hpp"
#include "utils/TtRecoFile.hpp"
#include "veto/Veto.hpp"
#include "veto/VetoTimeSaver.hpp"

class JRAF : public AlgBase {

public:

    JRAF(const std::string& name);

    ~JRAF() override = default;

    bool initialize() override;
    bool execute() override;
    bool finalize() override;

private:

    std::size_t m_iEvt = 0ul;
    TimeStamp m_tsEvt{0, 0};
    
    JM::NavBuffer* m_buf;
    RootInputSvc* m_iptSvc;

    // Reconstruction

    std::string m_recToolName;
	IRecMuonTool* m_recTool;
    std::string m_loaderName;
    std::string m_cdFillerName;
    std::string m_wpFillerName;
    std::string m_ttFillerName;
    Loader* m_loader;
    PmtTable m_pmtTable;
    Params m_params;

    // Muon selection variable

    double m_cd_muon_totq_thold = 10000.0;
    double m_cd_only_muon_totq_thold = 30000.0;
    double m_wp_muon_totq_thold = 600.0;
    TimeStamp m_cd_afterpulse_thold{0, 50000};
    TimeStamp m_wp_afterpulse_thold{0, 4000};
    TimeStamp m_cd_last_muon{0, 0};
    TimeStamp m_wp_last_muon{0, 0};
    
    TtRecoFile m_ttRecoFile;
    ContextFileTracker m_contextTracker;
    BeginningOfJobVetoTracker m_begOfJobVetoTrkr;
    MissingHeaderVetoTracker m_missHdrVetoTrkr;
    BigGapsVetoTracker m_bigGapsVetoTrkr;
    MuonVetoTracker m_muvetoTrkr;

    std::string m_ofilename;
    TFile* m_file;
    DAQTimeSaver m_daqTimeSaver;
    VetoTimeSaver m_vetoTimeSaver;
    std::vector<std::string> m_methods;
    std::vector<std::shared_ptr<Analysis>> m_analyses;

    TrackSaver m_trkSaver;
    FeatureSaver m_featureSaver;

	bool initBufSvc();
    bool initRecTool();
    bool initLoader();
    bool initAnalyses();

    calibration_context getCalibrationContext(const std::list<JM::CalibPmtChannel*>& clb_list);
    DetectorType getDetectorType(JM::EvtNavigator* nav);
    int getTtLayerId(double z);
    void addTtToTrack(std::vector<track>& tracks, const TimeStamp& curts);
    void addFeature(const std::vector<track>& tracks, const TimeStamp& curts, int run_id);

    void addTrack(RecTrks& rec_tracks, const std::string& method, double totq_cd, double totq_wp, const TimeStamp& ts, const track::loc& det, std::vector<track>& tracks);
    void addTrack(JM::CdTrackRecHeader* cdt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks);
    void addTrack(JM::WpRecHeader* wpt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks);
    void addTrack(JM::TtRecHeader* ttt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks);
    void addVertex(JM::OecHeader* oec_hdr, const std::string& method, const TimeStamp& ts, const calibration_context& calib, std::vector<vertex>& vertices);
    void addVertex(JM::CdVertexRecHeader* cdv_hdr, const std::string& method, const TimeStamp& ts, const calibration_context& calib, std::vector<vertex>& vertices);

};

#endif // JRAF_JRAF_HPP_