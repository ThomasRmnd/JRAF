#include "JRAF.hpp"

#include <cmath>
#include <numeric>

#include "SniperKernel/AlgFactory.h"
#include "SniperKernel/SniperLog.h"

#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdTrackRecHeader.h"
#include "Event/CdTriggerHeader.h"
#include "Event/OecHeader.h"
#include "Event/TtRecHeader.h"
#include "Event/WpCalibHeader.h"
#include "Event/WpRecHeader.h"
#include "Event/WpTriggerHeader.h"
#include "EvtNavigator/EvtNavHelper.h"

#include "analysis/AccidentalAnalysis.hpp"
#include "analysis/IBDAnalysis.hpp"
#include "analysis/MultiplicityAnalysis.hpp"
#include "analysis/NeutronAnalysis.hpp"
#include "event/EventCache.hpp"
#include "utils/NavBufferWrapper.hpp"

DECLARE_ALGORITHM(JRAF);

JRAF::JRAF(const std::string& name) : 
    AlgBase{name}
{

    declProp("EventBuilder", m_eventBuilderName = "EventBuilder");

    declProp("ContextPreviousFilename", m_contextTracker.prevctx);
    declProp("ContextNextFilename", m_contextTracker.nextctx);

    declProp("OutputFilename", m_ofilename = "output.root");
    declProp("RecoTrackOutputFilename", m_trkSaver.filename = "");
    declProp("FeatureOutputFilename", m_featureSaver.filename = "");
}

bool JRAF::initialize() {
    if (!initBufSvc()) return false;

    SniperPtr<RootInputSvc> iptSvc(getParent(), "InputSvc");
    if (iptSvc.invalid()) {
        LogError << "Can't find InputSvc\n";
        return false;
    }
    m_iptSvc = iptSvc.data();

    if (!m_trkSaver.init()) return false;
    if (!m_featureSaver.init()) return false;

    m_eventBuilder = tool<EventBuilder>(m_eventBuilderName);
    if (!m_eventBuilder) {
        LogError << "Failed to retrieve EventBuilder tool named " << m_eventBuilderName << '\n';
        return false;
    }
    if (!m_eventBuilder->initialize()) return false;

    if (!initAnalyses()) return false;

    LogInfo  << objName() << " initialized successfully\n"; 
    return true;
}

bool JRAF::initBufSvc() {
    SniperDataPtr<JM::NavBuffer> navBuf(getParent(), "/Event");
    if (navBuf.invalid()) {
        LogError << "Cannot get the NavBuffer at /Event\n";
        return false;
    }
    m_buf = navBuf.data();
    return true; 
}

bool JRAF::initAnalyses() {
    m_file = TFile::Open(m_ofilename.c_str(), "RECREATE");
    if (!m_file) {
        LogError << "Failed to create output file: " << m_ofilename << '\n';
        return false;
    }

    if (!m_daqTimeSaver.init()) return false;
    if (!m_vetoTimeSaver.init()) return false;

    m_methods = std::vector<std::string>{/* "Oec", */ "OMILREC", "MixedPhase", "OMILREC_JVtx" /* "JVertex" */};

    // m_analyses.push_back(std::make_shared<IBDAnalysis>("IBDAnalysis__Oec", "Oec"));
    m_analyses.push_back(std::make_shared<IBDAnalysis>("IBDAnalysis__OMILREC", "OMILREC"));
    m_analyses.push_back(std::make_shared<IBDAnalysis>("IBDAnalysis__MixedPhase", "MixedPhase"));
    m_analyses.push_back(std::make_shared<IBDAnalysis>("IBDAnalysis__OMILREC_JVtx", "OMILREC_JVtx"));
    m_analyses.push_back(std::make_shared<AccidentalAnalysis>("AccidentalAnalysis__OMILREC_JVtx", "OMILREC_JVtx"));
    m_analyses.push_back(std::make_shared<MultiplicityAnalysis>("MultiplicityAnalysis__OMILREC_JVtx", "OMILREC_JVtx"));
    m_analyses.push_back(std::make_shared<NeutronAnalysis>("NeutronAnalysis__OMILREC_JVtx", "OMILREC_JVtx"));

    for (std::shared_ptr<Analysis>& analysis : m_analyses) {
        if (!analysis->initialize()) return false;
    }
    return true;
}

// void JRAF::addFeature(const std::vector<track>& tracks, const TimeStamp& curts, int run_id) {
//     m_featureSaver.reset();
    
//     if (!m_ttRecoFile.find(curts)) return;
//     if (m_ttRecoFile.NTracks != 1) {
//         LogInfo << "Muon event is empty or a bundle considering TT (" << m_ttRecoFile.NTracks << " tracks, " << m_ttRecoFile.NTotPoints << " points)\n";
//         return;
//     }
//     if (m_ttRecoFile.NPoints[0] < 3) {
//         LogInfo << "Muon track has less than 3 points in the TT\n";
//         return;
//     }

//     std::unordered_set<int> layers_hit;
//     layers_hit.reserve(6);

//     for (int i = 0; i < m_ttRecoFile.NTotPoints; ++i) {
//         int lid = getTtLayerId(m_ttRecoFile.PointZ[i] + 26452.0);
//         if (lid >= 0) layers_hit.insert(lid);
//     }

//     if (layers_hit.size() < 3) {
//         LogInfo << "Muon track is not in three different layers of the TT\n";
//         return;
//     }

//     std::map<std::string, std::vector<std::vector<track>::const_iterator>> track_map;
//     track_map["CdWpTtChi2"] = {};
//     track_map["CdClassify"] = {};
//     track_map["WpBasic"] = {};
//     for (std::vector<track>::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
//         track_map[it->method].push_back(it);
//     }
//     if (track_map["CdWpTtChi2"].size() != 1 || track_map["CdClassify"].size() != 1 || track_map["WpBasic"].size() != 1) {
//         LogInfo << "Muon event is empty or a bundle considering CdWpTtChi2 or CdClassify or WpBasic\n";
//         return;
//     }
//     std::vector<track>::const_iterator trk_cdwpttchi2 = track_map["CdWpTtChi2"][0];
//     std::vector<track>::const_iterator trk_cdclassify = track_map["CdClassify"][0];
//     std::vector<track>::const_iterator trk_wpbasic = track_map["WpBasic"][0];

//     m_featureSaver.run_id = run_id;
//     m_featureSaver.sec = curts.GetSec();
//     m_featureSaver.nsec = curts.GetNanoSec();

//     m_featureSaver.iposx.push_back(trk_cdwpttchi2->ipos.x);
//     m_featureSaver.iposy.push_back(trk_cdwpttchi2->ipos.y);
//     m_featureSaver.iposz.push_back(trk_cdwpttchi2->ipos.z);
//     m_featureSaver.fposx.push_back(trk_cdwpttchi2->fpos.x);
//     m_featureSaver.fposy.push_back(trk_cdwpttchi2->fpos.y);
//     m_featureSaver.fposz.push_back(trk_cdwpttchi2->fpos.z);
//     m_featureSaver.chi2.push_back(trk_cdwpttchi2->quality);
//     m_featureSaver.det.push_back(0b001);

//     m_featureSaver.iposx.push_back(trk_wpbasic->ipos.x);
//     m_featureSaver.iposy.push_back(trk_wpbasic->ipos.y);
//     m_featureSaver.iposz.push_back(trk_wpbasic->ipos.z);
//     m_featureSaver.fposx.push_back(trk_wpbasic->fpos.x);
//     m_featureSaver.fposy.push_back(trk_wpbasic->fpos.y);
//     m_featureSaver.fposz.push_back(trk_wpbasic->fpos.z);
//     m_featureSaver.chi2.push_back(trk_wpbasic->quality);
//     m_featureSaver.det.push_back(0b010);

//     m_featureSaver.iposx.push_back(trk_cdclassify->ipos.x);
//     m_featureSaver.iposy.push_back(trk_cdclassify->ipos.y);
//     m_featureSaver.iposz.push_back(trk_cdclassify->ipos.z);
//     m_featureSaver.fposx.push_back(trk_cdclassify->fpos.x);
//     m_featureSaver.fposy.push_back(trk_cdclassify->fpos.y);
//     m_featureSaver.fposz.push_back(trk_cdclassify->fpos.z);
//     m_featureSaver.chi2.push_back(trk_cdclassify->quality);
//     m_featureSaver.det.push_back(0b011);

//     vec3 ipos(m_ttRecoFile.Coeff0[0], m_ttRecoFile.Coeff1[0], m_ttRecoFile.Coeff2[0] + 26452.0);
//     vec3 dir = unit(vec3{m_ttRecoFile.Coeff3[0], m_ttRecoFile.Coeff4[0], m_ttRecoFile.Coeff5[0]});
//     vec3 fpos = ipos - 2.0 * dot(ipos, dir) * dir;
//     m_featureSaver.iposx.push_back(ipos.x);
//     m_featureSaver.iposy.push_back(ipos.y);
//     m_featureSaver.iposz.push_back(ipos.z);
//     m_featureSaver.fposx.push_back(fpos.x);
//     m_featureSaver.fposy.push_back(fpos.y);
//     m_featureSaver.fposz.push_back(fpos.z);
//     m_featureSaver.chi2.push_back(m_ttRecoFile.Chi2[0]);
//     m_featureSaver.det.push_back(0b100);

//     for (const PmtProp& pmt : m_pmtTable) {
//         if (!pmt.used) continue;
//         if (pmt.type != Pmttype::_PMTINCH20 || (pmt.loc != 1 && pmt.loc != 2)) continue;
//         m_featureSaver.id.push_back(pmt.pmtid);
//         m_featureSaver.fht.push_back(pmt.fht);
//         m_featureSaver.totq.push_back(pmt.q);
//         m_featureSaver.q.push_back(static_cast<double>(pmt.hitq[0]));
//         m_featureSaver.nhit.push_back(pmt.hitq.size());
//     }

//     for (Int_t k = 0; k < m_ttRecoFile.NTotPoints; ++k) {
//         m_featureSaver.pointx.push_back(m_ttRecoFile.PointX[k]);
//         m_featureSaver.pointy.push_back(m_ttRecoFile.PointY[k]);
//         m_featureSaver.pointz.push_back(m_ttRecoFile.PointZ[k] + 26452.0);
//     }

//     m_featureSaver.fill();
// }

bool JRAF::execute() {
    LogInfo << "---------- Processing event by JRAF: " << ++m_iEvt << " ----------\n";

    JM::EvtNavigator* nav = m_buf->curEvt();
    if (!nav) {
        LogError << "EvtNavigator is nullptr\n";
        return false;
    }
    m_tsEvt = TimeStamp{nav->TimeStamp().GetTimeSpec()};
    int runId = nav->RunID();
    LogInfo << "TimeStamp: " << m_tsEvt << ", RunID: " << runId << '\n';

    if (!m_contextTracker.isTarget(m_iptSvc)) return true;

    // DEBUG --- Timing
    using clock = std::chrono::steady_clock;
    auto t_start = clock::now();
    // DEBUG --- Timing

    m_eventBuilder->build(m_buf);

    // DEBUG --- Timing
    auto t_after_load = clock::now();
    // DEBUG --- Timing

    m_daqTimeSaver.add(m_tsEvt, runId);

    if (m_begOfJobVetoTrkr.check(m_iEvt)) {
        if (!m_vetoTimeSaver.create(m_tsEvt, VetoType::BeginningOfJob, runId)) {
            LogError << "Failed to create veto time saver for beginning of job\n";
            return false;
        }
    }
    if (m_missHdrVetoTrkr.check(nav)) {
        if (!m_vetoTimeSaver.create(m_tsEvt, VetoType::MissingHeaders, runId)) {
            LogError << "Failed to create veto time saver for missing headers\n";
            return false;
        }
    }
    if (m_bigGapsVetoTrkr.check(nav)) {
        if (!m_vetoTimeSaver.create(m_tsEvt, VetoType::BigGaps, runId)) {
            LogError << "Failed to create veto time saver for big gaps\n";
            return false;
        }
    }
    if (m_muvetoTrkr.check(nav)) {
        if (!m_vetoTimeSaver.create_no_veto(m_tsEvt, VetoType::Muon, runId)) {
            LogError << "Failed to create veto time saver for muon\n";
            return false;
        }
    }

    m_trkSaver.fill(nav);
        
    if (m_vetoTimeSaver.inVeto(m_tsEvt)) return true;

    // DEBUG --- Timing
    auto t_before_context = clock::now();
    // DEBUG --- Timing
    
    EventContext events(m_buf, m_methods);

    // DEBUG --- Timing
    auto t_after_context = clock::now();
    // DEBUG --- Timing

    // --------------------
    // Per-analysis timing
    // --------------------
    std::vector<long long> per_analysis_ms;
    per_analysis_ms.reserve(m_analyses.size());
    for (std::shared_ptr<Analysis>& analysis : m_analyses) {
        auto t_before = clock::now();

        analysis->process(events.view(analysis->method()));

        auto t_after = clock::now();
        per_analysis_ms.push_back(
            std::chrono::duration_cast<std::chrono::milliseconds>(t_after - t_before).count()
        );
    }

    // DEBUG --- Timing
    auto t_after_analysis = clock::now();
    // DEBUG --- Timing

    auto t_load_ms  = std::chrono::duration_cast<std::chrono::milliseconds>(t_after_load - t_start).count();
    auto t_context_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_after_context - t_before_context).count();
    auto t_analysis_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_after_analysis - t_after_context).count();

    std::cout << "\n=== Timing report ===\n";
    std::cout << "1. Loading:  " << t_load_ms << " ms\n";
    std::cout << "2. Context: " << t_context_ms << " ms\n";
    std::cout << "3. Analysis: " << t_analysis_ms << " ms\n";
    std::cout << "=====================\n\n";

    std::cout << "-- Per-analysis breakdown --\n";
    for (size_t i = 0; i < m_analyses.size(); ++i) {
        std::cout << "   [" << m_analyses[i]->method() << "] "
                  << per_analysis_ms[i] << " ms\n";
    }
    std::cout << "=====================\n\n";

    return true;
}

bool JRAF::finalize() {
    if (m_eventBuilder && !m_eventBuilder->finalize()) return false;
    
    if (!m_trkSaver.save()) return false;
    if (!m_featureSaver.save()) return false;

    if (!m_file) return false;
    m_file->cd();
    if (!m_daqTimeSaver.write()) return false;
    if (!m_vetoTimeSaver.write()) return false;
    for (std::shared_ptr<Analysis>& ana : m_analyses) {
        if (!ana->write()) return false;
    }
    m_file->Close();

    LogInfo  << objName() << " finalized successfully\n"; 
    return true;
}