#include "AnalysisGroupC.hpp"

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

#include "analysis/CdWpCosmoStudy.hpp"
#include "analysis/FirstCrossCheckAnalysis.hpp"
#include "analysis/IBDWithCylindricalCut.hpp"
#include "analysis/IBDWithNeutronVetoStudy.hpp"
#include "analysis/MultiplicityWindowCut.hpp"
#include "analysis/NeutronVetoStudy.hpp"
#include "analysis/TtCosmoStudy.hpp"
#include "event/EventCache.hpp"

DECLARE_ALGORITHM(AnalysisGroupC);

AnalysisGroupC::AnalysisGroupC(const std::string& name) : 
    AlgBase{name}
{
    declProp("Loader", m_loaderName = "JointLoader");
    declProp("CdFiller", m_cdFillerName = "CdRangeFiller");
    declProp("WpFiller", m_wpFillerName = "WpRangeFiller");
    declProp("TtFiller", m_ttFillerName = "TtRangeFiller");

    declProp("RecTool", m_recToolName);
    declProp("ClassifyTool", m_classifyToolName);;

    declProp("TtRecoFilepath", m_ttRecoFile.filename = "");
    declProp("OutputFilename", m_ofilename = "output.root");
    declProp("TargetInputFilename", m_contextTracker.target);
}

bool AnalysisGroupC::initialize() {
    if (!initBufSvc()) return false;

    SniperPtr<RootInputSvc> iptSvc(getParent(), "InputSvc");
    if (iptSvc.invalid()) {
        LogError << "Can't find InputSvc." << std::endl;
        return false;
    }
    m_iptSvc = iptSvc.data();

    m_cd_last_muon = TimeStamp{0, 0};
    m_wp_last_muon = TimeStamp{0, 0};
    if (!m_ttRecoFile.init()) return false;
    if (!initLoader()) return false;
    if (!initRecTool()) return false;
    if (!initAnalyses()) return false;

    LogInfo  << objName() << " initialized successfully\n"; 
    return true;
}

bool AnalysisGroupC::initBufSvc() {
    SniperDataPtr<JM::NavBuffer> navBuf(getParent(), "/Event");
    if (navBuf.invalid()) {
        LogError << "Cannot get the NavBuffer @ /Event\n";
        return false;
    }
    m_buf = navBuf.data();
    return true; 
}

bool AnalysisGroupC::initLoader() {
    m_loader = tool<Loader>(m_loaderName);
    if (!m_loader) {
        LogError << "Failed to retrieve reconstruction tool named " << m_loaderName << '\n';
        return false;
    }
    RangeFiller<CdFillerTag>* cd_filler = tool<RangeFiller<CdFillerTag>>(m_cdFillerName);
    RangeFiller<WpFillerTag>* wp_filler = tool<RangeFiller<WpFillerTag>>(m_wpFillerName);
    RangeFiller<TtFillerTag>* tt_filler = tool<RangeFiller<TtFillerTag>>(m_ttFillerName);
    if (!m_loader->configure(&m_pmtTable, cd_filler, wp_filler, tt_filler)) return false;
	if (!m_loader->initialize()) return false;

    m_classifyLoader = tool<Loader>("BasicLoader");
    if (!m_classifyLoader) {
        LogError << "Failed to retrieve reconstruction tool named BasicLoader\n";
        return false;
    }
    wp_filler = tool<RangeFiller<WpFillerTag>>("AkiraWpRangeFiller");
    if (!m_classifyLoader->configure(&m_classifyPmtTable, nullptr, wp_filler, nullptr)) return false;
    if (!m_classifyLoader->initialize()) return false;

    return true;
}

bool AnalysisGroupC::initRecTool() {
    m_recTool = tool<IRecMuonTool>(m_recToolName);
    if (!m_recTool) {
        LogError << "Failed to retrieve reconstruction tool named " << m_recToolName << '\n';
        return false;
    }
    if (!m_recTool->configure(&m_params, &m_pmtTable)) return false;
	// if (!dynamic_cast<ToolBase*>(m_recTool)->initialize()) return false;

    m_classifyTool = tool<IRecMuonTool>(m_classifyToolName);
    if (!m_classifyTool) {
        LogError << "Failed to retrieve classification tool named " << m_classifyToolName << '\n';
        return false;
    }
    if (!m_classifyTool->configure(&m_classifyParams, &m_classifyPmtTable)) return false;
    // if (!dynamic_cast<ToolBase*>(m_classifyTool)->initialize()) return false;
    return true;
}

bool AnalysisGroupC::initAnalyses() {
    m_file = TFile::Open(m_ofilename.c_str(), "RECREATE");
    if (!m_file) {
        LogError << "Failed to create output file: " << m_ofilename << '\n';
        return false;
    }

    m_daq_tree = new TTree("DAQTree", "DAQTree");
    if (!m_daq_tree) {
        LogError << "Failed to create TTree DAQTree\n";
        return false;
    }
    m_daq_tree->Branch("daq_sec", &m_daq_sec);
    m_daq_tree->Branch("daq_nsec", &m_daq_nsec);
    m_daq_tree->Branch("muveto_sec", &m_muveto_sec);
    m_daq_tree->Branch("muveto_nsec", &m_muveto_nsec);
    m_daq_sec = 0l;
    m_daq_nsec = 0;
    m_muveto_sec = 0l;
    m_muveto_nsec = 0;

    // m_analyses.push_back(std::make_shared<FirstCrossCheckAnalysis>("FirstCrossCheckAnalysis__Oec", "Oec"));
    // m_analyses.push_back(std::make_shared<IBDWithCylindricalCut>("IBDWithCylindricalCut_3m__Oec", "Oec", 3000.0));
    // m_analyses.push_back(std::make_shared<IBDWithCylindricalCut>("IBDWithCylindricalCut_5m__Oec", "Oec", 5000.0));
    // m_analyses.push_back(std::make_shared<MultiplicityWindowCut>("MultiplicityWindowCut__Oec", "Oec"));
    // m_analyses.push_back(std::make_shared<TtCosmoStudy>("TtCosmoStudy_sig__Oec", "Oec", TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}));
    // m_analyses.push_back(std::make_shared<TtCosmoStudy>("TtCosmoStudy_bkg__Oec", "Oec", TimeStamp{0, -1200000000}, TimeStamp{0, -5000000}));
    // m_analyses.push_back(std::make_shared<CdWpCosmoStudy>("CdWpCosmoStudy_All__Oec", "Oec", 40000.0, TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}));
    // m_analyses.push_back(std::make_shared<CdWpCosmoStudy>("CdWpCosmoStudy_sig__Oec", "Oec", 3000.0, TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}));
    // m_analyses.push_back(std::make_shared<CdWpCosmoStudy>("CdWpCosmoStudy_bkg__Oec", "Oec", 3000.0, TimeStamp{0, -1200000000}, TimeStamp{0, -5000000}));
    // m_analyses.push_back(std::make_shared<IBDWithNeutronVetoStudy>("IBDWithNeutronVetoStudy__Oec", "Oec"));
    // m_analyses.push_back(std::make_shared<NeutronVetoStudy>("NeutronVetoStudy_3m_1_5s__Oec", "Oec", 3000.0, TimeStamp{0, 1500000000}));
    // m_analyses.push_back(std::make_shared<NeutronVetoStudy>("NeutronVetoStudy_All__Oec", "Oec", 40000.0, TimeStamp{0, 2000000000}));

    m_analyses.push_back(std::make_shared<FirstCrossCheckAnalysis>("FirstCrossCheckAnalysis__OMILREC", "OMILREC"));
    m_analyses.push_back(std::make_shared<IBDWithCylindricalCut>("IBDWithCylindricalCut_3m__OMILREC", "OMILREC", 3000.0));
    m_analyses.push_back(std::make_shared<IBDWithCylindricalCut>("IBDWithCylindricalCut_5m__OMILREC", "OMILREC", 5000.0));
    m_analyses.push_back(std::make_shared<MultiplicityWindowCut>("MultiplicityWindowCut__OMILREC", "OMILREC"));
    m_analyses.push_back(std::make_shared<TtCosmoStudy>("TtCosmoStudy_sig__OMILREC", "OMILREC", TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}));
    m_analyses.push_back(std::make_shared<TtCosmoStudy>("TtCosmoStudy_bkg__OMILREC", "OMILREC", TimeStamp{0, -1200000000}, TimeStamp{0, -5000000}));
    m_analyses.push_back(std::make_shared<CdWpCosmoStudy>("CdWpCosmoStudy_All__OMILREC", "OMILREC", 40000.0, TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}));
    m_analyses.push_back(std::make_shared<CdWpCosmoStudy>("CdWpCosmoStudy_sig__OMILREC", "OMILREC", 3000.0, TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}));
    m_analyses.push_back(std::make_shared<CdWpCosmoStudy>("CdWpCosmoStudy_bkg__OMILREC", "OMILREC", 3000.0, TimeStamp{0, -1200000000}, TimeStamp{0, -5000000}));
    m_analyses.push_back(std::make_shared<IBDWithNeutronVetoStudy>("IBDWithNeutronVetoStudy__OMILREC", "OMILREC"));
    m_analyses.push_back(std::make_shared<NeutronVetoStudy>("NeutronVetoStudy_3m_1_5s__OMILREC", "OMILREC", 3000.0, TimeStamp{0, 1500000000}));
    m_analyses.push_back(std::make_shared<NeutronVetoStudy>("NeutronVetoStudy_All__OMILREC", "OMILREC", 40000.0, TimeStamp{0, 2000000000}));

    m_analyses.push_back(std::make_shared<FirstCrossCheckAnalysis>("FirstCrossCheckAnalysis__MixedPhase", "MixedPhase"));
    m_analyses.push_back(std::make_shared<IBDWithCylindricalCut>("IBDWithCylindricalCut_3m__MixedPhase", "MixedPhase", 3000.0));
    m_analyses.push_back(std::make_shared<IBDWithCylindricalCut>("IBDWithCylindricalCut_5m__MixedPhase", "MixedPhase", 5000.0));
    m_analyses.push_back(std::make_shared<MultiplicityWindowCut>("MultiplicityWindowCut__MixedPhase", "MixedPhase"));
    m_analyses.push_back(std::make_shared<TtCosmoStudy>("TtCosmoStudy_sig__MixedPhase", "MixedPhase", TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}));
    m_analyses.push_back(std::make_shared<TtCosmoStudy>("TtCosmoStudy_bkg__MixedPhase", "MixedPhase", TimeStamp{0, -1200000000}, TimeStamp{0, -5000000}));
    m_analyses.push_back(std::make_shared<CdWpCosmoStudy>("CdWpCosmoStudy_All__MixedPhase", "MixedPhase", 40000.0, TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}));
    m_analyses.push_back(std::make_shared<CdWpCosmoStudy>("CdWpCosmoStudy_sig__MixedPhase", "MixedPhase", 3000.0, TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}));
    m_analyses.push_back(std::make_shared<CdWpCosmoStudy>("CdWpCosmoStudy_bkg__MixedPhase", "MixedPhase", 3000.0, TimeStamp{0, -1200000000}, TimeStamp{0, -5000000}));
    m_analyses.push_back(std::make_shared<IBDWithNeutronVetoStudy>("IBDWithNeutronVetoStudy__MixedPhase", "MixedPhase"));
    m_analyses.push_back(std::make_shared<NeutronVetoStudy>("NeutronVetoStudy_3m_1_5s__MixedPhase", "MixedPhase", 3000.0, TimeStamp{0, 1500000000}));
    m_analyses.push_back(std::make_shared<NeutronVetoStudy>("NeutronVetoStudy_All__MixedPhase", "MixedPhase", 40000.0, TimeStamp{0, 2000000000}));

    m_analyses.push_back(std::make_shared<FirstCrossCheckAnalysis>("FirstCrossCheckAnalysis__JVertex", "JVertex"));
    m_analyses.push_back(std::make_shared<IBDWithCylindricalCut>("IBDWithCylindricalCut_3m__JVertex", "JVertex", 3000.0));
    m_analyses.push_back(std::make_shared<IBDWithCylindricalCut>("IBDWithCylindricalCut_5m__JVertex", "JVertex", 5000.0));
    m_analyses.push_back(std::make_shared<MultiplicityWindowCut>("MultiplicityWindowCut__JVertex", "JVertex"));
    m_analyses.push_back(std::make_shared<TtCosmoStudy>("TtCosmoStudy_sig__JVertex", "JVertex", TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}));
    m_analyses.push_back(std::make_shared<TtCosmoStudy>("TtCosmoStudy_bkg__JVertex", "JVertex", TimeStamp{0, -1200000000}, TimeStamp{0, -5000000}));
    m_analyses.push_back(std::make_shared<CdWpCosmoStudy>("CdWpCosmoStudy_All__JVertex", "JVertex", 40000.0, TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}));
    m_analyses.push_back(std::make_shared<CdWpCosmoStudy>("CdWpCosmoStudy_sig__JVertex", "JVertex", 3000.0, TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}));
    m_analyses.push_back(std::make_shared<CdWpCosmoStudy>("CdWpCosmoStudy_bkg__JVertex", "JVertex", 3000.0, TimeStamp{0, -1200000000}, TimeStamp{0, -5000000}));
    m_analyses.push_back(std::make_shared<IBDWithNeutronVetoStudy>("IBDWithNeutronVetoStudy__JVertex", "JVertex"));
    m_analyses.push_back(std::make_shared<NeutronVetoStudy>("NeutronVetoStudy_3m_1_5s__JVertex", "JVertex", 3000.0, TimeStamp{0, 1500000000}));
    m_analyses.push_back(std::make_shared<NeutronVetoStudy>("NeutronVetoStudy_All__JVertex", "JVertex", 40000.0, TimeStamp{0, 2000000000}));

    for (std::shared_ptr<Analysis>& analysis : m_analyses) {
        if (!analysis->initialize()) return false;
    }
    return true;
}

void AnalysisGroupC::addTrack(RecTrks& rec_tracks, const std::string& method, const TimeStamp& ts, const track::loc& det, std::vector<track>& tracks) {
    for (int k = 0; k < rec_tracks.size(); ++k) {
        tracks.push_back(track{
            method, vec3{rec_tracks.getStart(k)}, vec3{rec_tracks.getEnd(k)}, rec_tracks.getNPE(k), ts, det, rec_tracks.getQuality(k)
        });
    }
}

void AnalysisGroupC::addTrack(JM::CdTrackRecHeader* cdt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks) {
    if (!cdt_hdr || !cdt_hdr->event()) return;
    const std::vector<JM::RecTrack*>& rec_tracks = cdt_hdr->event()->tracks();
    for (JM::RecTrack* t : rec_tracks) {
        tracks.push_back(track{
            method, *t, ts, track::loc::cd
        });
    }
}

void AnalysisGroupC::addTrack(JM::WpRecHeader* wpt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks) {
    if (!wpt_hdr || !wpt_hdr->event()) return;
    const std::vector<JM::RecTrack*>& rec_tracks = wpt_hdr->event()->tracks();
    for (JM::RecTrack* t : rec_tracks) {
        tracks.push_back(track{
            method, *t, ts, track::loc::cd
        });
    }
}

void AnalysisGroupC::addTrack(JM::TtRecHeader* ttt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks) {
    if (!ttt_hdr || !ttt_hdr->event()) return;
    JM::TtRecEvt* ttt_evt = ttt_hdr->event();
    for (int k = 0; k < ttt_evt->nTracks(); ++k) {
        vec3 ipos{ttt_evt->Coeff0()[k], ttt_evt->Coeff1()[k], ttt_evt->Coeff2()[k]};
        vec3 dir = unit(vec3{ttt_evt->Coeff3()[k], ttt_evt->Coeff4()[k], ttt_evt->Coeff5()[k]});
        vec3 fpos = ipos - 2.0 * dot(ipos, dir) * dir;
        tracks.push_back(track{
            method, ipos, fpos, 0.0, ts, track::loc::tt, ttt_evt->Chi2()[k]
        });
    }
}

void AnalysisGroupC::addVertex(JM::OecHeader* oec_hdr, const std::string& method, const TimeStamp& ts, double totq, std::vector<vertex>& vertices) {
    if (!oec_hdr || !oec_hdr->event("JM::OecEvt")) return;
    JM::OecEvt* oec_evt = dynamic_cast<JM::OecEvt*>(oec_hdr->event("JM::OecEvt"));
    vertices.push_back(vertex{
        "Oec", vec3{oec_evt->getVertexX(), oec_evt->getVertexY(), oec_evt->getVertexZ()}, oec_evt->getEnergy(), totq, ts, "Unknown"
    });
}

void AnalysisGroupC::addVertex(JM::CdVertexRecHeader* cdv_hdr, const std::string& method, const TimeStamp& ts, double totq, std::vector<vertex>& vertices) {
    if (!cdv_hdr || !cdv_hdr->event()) return;
    const std::vector<JM::RecVertex*>& rec_vertices = cdv_hdr->event()->vertices();
    for (JM::RecVertex* v : rec_vertices) {
        vertices.push_back(vertex{
            method, vec3{v->x(), v->y(), v->z()}, v->energy(), totq, ts, "Unknown"
        });
    }
}

bool AnalysisGroupC::execute() {
    LogInfo << "---------- Processing event by AnalysisGroupC: " << ++m_iEvt << " ----------\n";

    JM::EvtNavigator* nav = m_buf->curEvt();
    if (!nav) {
        LogError << "EvtNavigator is nullptr\n";
        return false;
    }
    m_tsEvt = TimeStamp{nav->TimeStamp().GetTimeSpec()};

    LogInfo << "TimeStamp: " << m_tsEvt << '\n';

    if (!m_contextTracker.isTarget(m_iptSvc)) return true;
    if (m_iEvt == 1ul) {
        m_targetIsFirst = true;
        m_targetFirstTs = m_tsEvt;
    }
    LogInfo << "Target input is first " << m_targetIsFirst << '\n';
    if (m_targetIsFirst) {
        TimeStamp ts_diff = m_tsEvt - m_targetFirstTs;
        if (ts_diff <= TimeStamp{0, 5000000}) return true;
    }

    NavBufferWrapper bufwrap(*m_buf);
    for (; bufwrap.current() != bufwrap.end(); bufwrap.next()) {
        if (EventCache::contains(bufwrap.curEvt())) continue;

        std::shared_ptr<Event> evt = std::make_shared<Event>();

        if (!m_loader->load(&bufwrap)) return false;
        if (!m_classifyLoader->load(&bufwrap)) return false;

        TimeStamp curts{bufwrap.curEvt()->TimeStamp().GetTimeSpec()};

        double totq_cd = std::accumulate(m_pmtTable.begin(), m_pmtTable.end(), 0.0, 
            [](double sum, const PmtProp& pmt) { return sum + ( (pmt.used && (pmt.type & PmtType::PMT_20INCH) == pmt.type) ? pmt.q : 0.0 ); } 
        );
        double totq_wp = std::accumulate(m_pmtTable.begin(), m_pmtTable.end(), 0.0, 
            [](double sum, const PmtProp& pmt) { return sum + ( (pmt.used && (pmt.type & PmtType::PMT_WP) == pmt.type) ? pmt.q : 0.0 ); } 
        );

        LogInfo << "TotQ: CD = " << totq_cd << ", WP = " << totq_wp << '\n';

        bool is_possibly_cd_muon = false;
        bool is_possibly_wp_muon = false;

        if (
            totq_cd >= m_cd_muon_totq_thold && 
            totq_wp >= m_wp_muon_totq_thold && 
            curts - m_cd_last_muon > m_cd_afterpulse_thold &&
            curts - m_wp_last_muon > m_wp_afterpulse_thold
        ) {
            m_cd_last_muon = curts;
            m_wp_last_muon = curts;
            is_possibly_cd_muon = true;
            is_possibly_wp_muon = true;
        }
        else if (
            totq_cd < m_cd_muon_totq_thold && 
            totq_wp >= m_wp_muon_totq_thold &&
            curts - m_wp_last_muon > m_wp_afterpulse_thold
        ) {
            m_wp_last_muon = curts;
            is_possibly_wp_muon = true;
        }
        else if (
            totq_cd >= m_cd_muon_totq_thold && 
            totq_wp < m_wp_muon_totq_thold &&
            (curts - m_cd_last_muon > TimeStamp{0, 2000000} || curts - m_wp_last_muon > TimeStamp{0, 2000000})
        ) {
            m_cd_last_muon = curts;
            is_possibly_cd_muon = true;
        }

        LogInfo << "Is possibly CD muon: " << is_possibly_cd_muon << ", is possibly WP muon: " << is_possibly_wp_muon << '\n';

        std::vector<track> tracks;
        if (is_possibly_cd_muon) {
            JM::CdTrackRecHeader* basic_cdt_hdr = JM::getHeaderObject<JM::CdTrackRecHeader>(bufwrap.curEvt());
            addTrack(basic_cdt_hdr, "CdBasic", curts, tracks);
            JM::CdTrackRecHeader* classify_cdt_hdr = JM::getHeaderObject<JM::CdTrackRecHeader>(bufwrap.curEvt(), "/Event/CdTrackRecClassify");
            addTrack(classify_cdt_hdr, "CdClassify", curts, tracks);

            RecTrks rtrks;
            if (!m_recTool->reconstruct(&rtrks)) {
                LogWarn << "Could not reconstruct the event with reconstruction tool\n";
            }
            addTrack(rtrks, "CdWpTtChi2", curts, track::loc::cd, tracks);
        }
        else if (is_possibly_wp_muon) {
            JM::WpRecHeader* basic_wpt_hdr = JM::getHeaderObject<JM::WpRecHeader>(bufwrap.curEvt());
            addTrack(basic_wpt_hdr, "WpBasic", curts, tracks);
            JM::WpRecHeader* classify_wpt_hdr = JM::getHeaderObject<JM::WpRecHeader>(bufwrap.curEvt(), "/Event/WpTrackRecClassify");
            addTrack(classify_wpt_hdr, "CdClassify", curts, tracks);
            
            RecTrks rtrks;
            if (!m_classifyTool->reconstruct(&rtrks)) {
                LogWarn << "Could not classify the event with classification tool\n";
            }
            addTrack(rtrks, "WpClassify", curts, track::loc::wp, tracks);
        }
        if ( (is_possibly_cd_muon || is_possibly_wp_muon) && tracks.empty() ) {
            tracks.push_back(track{"CdBasic", vec3{0.0, 0.0, 20000.0}, vec3{0.0, 0.0, -20000.0}, 0.0, curts, track::loc::cd, -1.0});
        }

        if (m_ttRecoFile.find(curts)) {
            if (m_ttRecoFile.ntracks > 100) {
                LogWarn << "More than 100 tracks reconstructed by the TT!\n";
            }

            bool is_good_tt_reco = false;
            if (m_ttRecoFile.ntracks >= 1) {
                for (Int_t k = 0; k < m_ttRecoFile.ntracks; ++k) {
                    if (m_ttRecoFile.npts[k] >= 3) {
                        is_good_tt_reco = true;
                        break;
                    }
                }
            }

            if (is_good_tt_reco) {
                for (Int_t k = 0; k < std::min(m_ttRecoFile.ntracks, 100); ++k) {
                    if (m_ttRecoFile.npts[k] < 3) continue;
                    vec3 ipos{m_ttRecoFile.coeff0[k], m_ttRecoFile.coeff1[k], m_ttRecoFile.coeff2[k]};
                    vec3 dir = unit(vec3{m_ttRecoFile.coeff3[k], m_ttRecoFile.coeff4[k], m_ttRecoFile.coeff5[k]});
                    vec3 fpos = ipos - 2.0 * dot(ipos, dir) * dir;
                    tracks.push_back(track{
                        "Tt", ipos, fpos, 0.0, curts, track::loc::tt, m_ttRecoFile.chi2[k]
                    });
                }
            }
        }

        std::vector<vertex> vertices;
        JM::CdVertexRecHeader* basic_cdv_hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(nav);
        addVertex(basic_cdv_hdr, "Basic", curts, totq_cd, vertices);
        JM::CdVertexRecHeader* jvertex_cdv_hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(nav, "/Event/CdVertexRecJVertex");
        addVertex(jvertex_cdv_hdr, "JVertex", curts, totq_cd, vertices);
        JM::CdVertexRecHeader* mixedphase_cdv_hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(nav, "/Event/CdVertexRecMixedPhase");
        addVertex(mixedphase_cdv_hdr, "MixedPhase", curts, totq_cd, vertices);
        JM::CdVertexRecHeader* omilrec_cdv_hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(nav, "/Event/CdVertexRecOMILREC");
        addVertex(omilrec_cdv_hdr, "OMILREC", curts, totq_cd, vertices);

        evt->tracks = tracks;
        evt->vertices = vertices;
        EventCache::insert(curts, evt);
        LogInfo << *evt << '\n';
    }

    JM::OecHeader* oec_hdr = JM::getHeaderObject<JM::OecHeader>(nav);
    JM::CdLpmtCalibHeader* cd_lpmt_calib_hdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
    // JM::CdTriggerHeader* cd_trig_hdr = JM::getHeaderObject<JM::CdTriggerHeader>(nav);
    JM::WpCalibHeader* wp_calib_hdr = JM::getHeaderObject<JM::WpCalibHeader>(nav);
    // JM::WpTriggerHeader* wp_trig_hdr = JM::getHeaderObject<JM::WpTriggerHeader>(nav);
        
    if (!oec_hdr || (!cd_lpmt_calib_hdr && !wp_calib_hdr) /* || (!cd_trig_hdr && !wp_trig_hdr) */) {
        m_vetoTs = m_tsEvt;
    }
    TimeStamp ts_diff = m_tsEvt - m_vetoTs;
    if (ts_diff <= TimeStamp{0, 5000000}) return true;

    for (std::shared_ptr<Analysis>& analysis : m_analyses) {
        analysis->process(m_buf);
    }
    if (m_buf->begin() <= m_buf->current() - 1l) {
        JM::EvtNavigator* prv_nav = (m_buf->current() - 1l)->get();
        TimeStamp prv_ts{prv_nav->TimeStamp().GetTimeSpec()};
        TimeStamp daq_ts{m_daq_sec, m_daq_nsec};
        daq_ts.Add(m_tsEvt - prv_ts);
        m_daq_sec = daq_ts.GetSec();
        m_daq_nsec = daq_ts.GetNanoSec();
    }
    std::shared_ptr<Event> evt = EventCache::load(nav);
    if (!evt->tracks.empty()) {
        TimeStamp muveto_ts{m_muveto_sec, m_muveto_nsec};
        muveto_ts.Add(TimeStamp{0, 5000000});
        m_muveto_sec = muveto_ts.GetSec();
        m_muveto_nsec = muveto_ts.GetNanoSec();
    }

    return true;
}

bool AnalysisGroupC::finalize() {
    if (m_loader && !m_loader->finalize()) return false;
    if (m_classifyLoader && !m_classifyLoader->finalize()) return false;
    if (m_recTool && !(dynamic_cast<ToolBase*>(m_recTool))->finalize()) return false;
    
    if (!m_file) return false;
    m_daq_tree->Fill();
    m_file->cd();
    m_daq_tree->Write();
    for (std::shared_ptr<Analysis>& ana : m_analyses) {
        if (!ana->write()) return false;
    }
    m_file->Close();

    LogInfo  << objName() << " finalized successfully\n"; 
    return true;
}