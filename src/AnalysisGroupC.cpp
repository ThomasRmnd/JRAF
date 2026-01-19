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

#include "analysis/IBDAnalysis.hpp"
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

    declProp("TtRecoFilepath", m_ttRecoFile.filename = "");
    declProp("ContextPreviousFilename", m_contextTracker.prevctx);
    declProp("ContextNextFilename", m_contextTracker.nextctx);

    declProp("OutputFilename", m_ofilename = "output.root");
    declProp("RecoTrackOutputFilename", m_trkSaver.filename = "");
    declProp("FeatureOutputFilename", m_featureSaver.filename = "");
}

bool AnalysisGroupC::initialize() {
    if (!initBufSvc()) return false;

    SniperPtr<RootInputSvc> iptSvc(getParent(), "InputSvc");
    if (iptSvc.invalid()) {
        LogError << "Can't find InputSvc." << std::endl;
        return false;
    }
    m_iptSvc = iptSvc.data();

    if (!m_ttRecoFile.init()) return false;
    if (!m_trkSaver.init()) return false;
    if (!m_featureSaver.init()) return false;
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
    return true;
}

bool AnalysisGroupC::initAnalyses() {
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
            method, *t, ts, track::loc::wp
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

void AnalysisGroupC::addVertex(JM::OecHeader* oec_hdr, const std::string& method, const TimeStamp& ts, const calibration_context& calib, std::vector<vertex>& vertices) {
    if (!oec_hdr || !oec_hdr->event("JM::OecEvt")) return;
    JM::OecEvt* oec_evt = dynamic_cast<JM::OecEvt*>(oec_hdr->event("JM::OecEvt"));
    vertices.push_back(vertex{
        method, vec3{oec_evt->getVertexX(), oec_evt->getVertexY(), oec_evt->getVertexZ()}, oec_evt->getEnergy(), ts, calib, "Unknown"
    });
}

void AnalysisGroupC::addVertex(JM::CdVertexRecHeader* cdv_hdr, const std::string& method, const TimeStamp& ts, const calibration_context& calib, std::vector<vertex>& vertices) {
    if (!cdv_hdr || !cdv_hdr->event()) return;
    const std::vector<JM::RecVertex*>& rec_vertices = cdv_hdr->event()->vertices();
    for (JM::RecVertex* v : rec_vertices) {
        vertices.push_back(vertex{
            method, vec3{v->x(), v->y(), v->z()}, v->energy(), ts, calib, "Unknown"
        });
    }
}

int AnalysisGroupC::getTtLayerId(double z) {
    if (24000.0 <= z && z <= 25000.0) return 0;  // main
    if (25500.0 <= z && z <= 26500.0) return 1;  // main
    if (27000.0 <= z && z <= 28000.0) return 2;  // main
    if (30000.0 <= z && z <= 30200.0) return 3;  // chimney
    if (30200.0 <= z && z <= 30400.0) return 4;  // chimney
    if (30400.0 <= z && z <= 30600.0) return 5;  // chimney

    return -1; // not inside any valid layer
};

void AnalysisGroupC::addTtToTrack(std::vector<track>& tracks, const TimeStamp& curts) {
    if (!m_ttRecoFile.find(curts)) return;
    
    if (m_ttRecoFile.NTracks < 1) {
        LogInfo << "No TT events found\n";
        return;
    }
    if (m_ttRecoFile.NTracks > 20) {
        LogWarn << "More than 20 tracks reconstructed by the TT!\n";
    }

    for (Int_t k = 0; k < std::min(m_ttRecoFile.NTracks, 20); ++k) {
        if (m_ttRecoFile.NPoints[k] < 3) continue;
        vec3 ipos{m_ttRecoFile.Coeff0[k], m_ttRecoFile.Coeff1[k], m_ttRecoFile.Coeff2[k] + 26452.0};
        vec3 dir = unit(vec3{m_ttRecoFile.Coeff3[k], m_ttRecoFile.Coeff4[k], m_ttRecoFile.Coeff5[k]});
        vec3 fpos = ipos - 2.0 * dot(ipos, dir) * dir;
        tracks.push_back(track{
            "Tt", ipos, fpos, 0.0, curts, track::loc::tt, m_ttRecoFile.Chi2[k]
        });
    }
}

void AnalysisGroupC::addFeature(const std::vector<track>& tracks, const TimeStamp& curts, int run_id) {
    m_featureSaver.reset();
    
    if (!m_ttRecoFile.find(curts)) return;
    if (m_ttRecoFile.NTracks != 1) {
        LogInfo << "Muon event is empty or a bundle considering TT (" << m_ttRecoFile.NTracks << " tracks, " << m_ttRecoFile.NTotPoints << " points)\n";
        return;
    }
    if (m_ttRecoFile.NPoints[0] < 3) {
        LogInfo << "Muon track has less than 3 points in the TT\n";
        return;
    }

    std::unordered_set<int> layers_hit;
    layers_hit.reserve(6);

    for (int i = 0; i < m_ttRecoFile.NTotPoints; ++i) {
        int lid = getTtLayerId(m_ttRecoFile.PointZ[i] + 26452.0);
        if (lid >= 0) layers_hit.insert(lid);
    }

    if (layers_hit.size() < 3) {
        LogInfo << "Muon track is not in three different layers of the TT\n";
        return;
    }

    std::map<std::string, std::vector<std::vector<track>::const_iterator>> track_map;
    track_map["CdWpTtChi2"] = {};
    track_map["CdClassify"] = {};
    for (std::vector<track>::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
        track_map[it->method].push_back(it);
    }
    if (track_map["CdWpTtChi2"].size() != 1 || track_map["CdClassify"].size() != 1) {
        LogInfo << "Muon event is empty or a bundle considering CdWpTtChi2 or CdClassify\n";
        return;
    }
    std::vector<track>::const_iterator trk_cdwpttchi2 = track_map["CdWpTtChi2"][0];
    std::vector<track>::const_iterator trk_cdclassify = track_map["CdClassify"][0];

    m_featureSaver.run_id = run_id;
    m_featureSaver.sec = curts.GetSec();
    m_featureSaver.nsec = curts.GetNanoSec();

    m_featureSaver.iposx.push_back(trk_cdwpttchi2->ipos.x);
    m_featureSaver.iposy.push_back(trk_cdwpttchi2->ipos.y);
    m_featureSaver.iposz.push_back(trk_cdwpttchi2->ipos.z);
    m_featureSaver.fposx.push_back(trk_cdwpttchi2->fpos.x);
    m_featureSaver.fposy.push_back(trk_cdwpttchi2->fpos.y);
    m_featureSaver.fposz.push_back(trk_cdwpttchi2->fpos.z);
    m_featureSaver.chi2.push_back(trk_cdwpttchi2->quality);
    m_featureSaver.det.push_back(1 << 0);

    m_featureSaver.iposx.push_back(trk_cdclassify->ipos.x);
    m_featureSaver.iposy.push_back(trk_cdclassify->ipos.y);
    m_featureSaver.iposz.push_back(trk_cdclassify->ipos.z);
    m_featureSaver.fposx.push_back(trk_cdclassify->fpos.x);
    m_featureSaver.fposy.push_back(trk_cdclassify->fpos.y);
    m_featureSaver.fposz.push_back(trk_cdclassify->fpos.z);
    m_featureSaver.chi2.push_back(trk_cdclassify->quality);
    m_featureSaver.det.push_back(1 << 1);

    vec3 ipos(m_ttRecoFile.Coeff0[0], m_ttRecoFile.Coeff1[0], m_ttRecoFile.Coeff2[0] + 26452.0);
    vec3 dir = unit(vec3{m_ttRecoFile.Coeff3[0], m_ttRecoFile.Coeff4[0], m_ttRecoFile.Coeff5[0]});
    vec3 fpos = ipos - 2.0 * dot(ipos, dir) * dir;
    m_featureSaver.iposx.push_back(ipos.x);
    m_featureSaver.iposy.push_back(ipos.y);
    m_featureSaver.iposz.push_back(ipos.z);
    m_featureSaver.fposx.push_back(fpos.x);
    m_featureSaver.fposy.push_back(fpos.y);
    m_featureSaver.fposz.push_back(fpos.z);
    m_featureSaver.chi2.push_back(m_ttRecoFile.Chi2[0]);
    m_featureSaver.det.push_back(1 << 2);

    for (const PmtProp& pmt : m_pmtTable) {
        if (!pmt.used) continue;
        if ( ( (pmt.type & PmtType::PMT_20INCH) != pmt.type ) && ( (pmt.type & PmtType::PMT_WP) != pmt.type ) ) continue;
        m_featureSaver.id.push_back(pmt.pmtid);
        m_featureSaver.fht.push_back(pmt.fht);
        m_featureSaver.totq.push_back(pmt.q);
        m_featureSaver.q.push_back(static_cast<double>(pmt.hitq[0]));
        m_featureSaver.nhit.push_back(pmt.hitq.size());
    }

    for (Int_t k = 0; k < m_ttRecoFile.NTotPoints; ++k) {
        m_featureSaver.pointx.push_back(m_ttRecoFile.PointX[k]);
        m_featureSaver.pointy.push_back(m_ttRecoFile.PointY[k]);
        m_featureSaver.pointz.push_back(m_ttRecoFile.PointZ[k] + 26452.0);
    }

    m_featureSaver.fill();
}

bool AnalysisGroupC::execute() {
    LogInfo << "---------- Processing event by AnalysisGroupC: " << ++m_iEvt << " ----------\n";

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

    NavBufferWrapper bufwrap(*m_buf);
    for (; bufwrap.current() != bufwrap.end(); bufwrap.next()) {
        if (EventCache::contains(bufwrap.curEvt())) continue;

        std::shared_ptr<Event> evt = std::make_shared<Event>();

        if (!m_loader->load(&bufwrap)) return false;

        TimeStamp curts{bufwrap.curEvt()->TimeStamp().GetTimeSpec()};

        calibration_context calib;
        double totq_wp = 0.0;
        for (PmtTable::const_iterator it = m_pmtTable.begin(); it != m_pmtTable.end(); ++it) {
            if (!it->used) continue;
            if ( (it->type & PmtType::PMT_20INCH) == it->type ) {
                calib.totq += it->q;
                calib.meant += it->fht;
                ++calib.npmt;
                calib.nhit += it->hittime.size();
                if (it->q < calib.minq) calib.minq = it->q;
                if (it->q > calib.maxq) calib.maxq = it->q;
            }
            else if ( (it->type & PmtType::PMT_WP) == it->type ) {
                totq_wp += it->q;
            }
        }
        if (calib.npmt > 0) {
            calib.meanq = calib.totq / calib.npmt;
            calib.meant = calib.meant / calib.npmt;
            calib.meanhit = static_cast<double>(calib.nhit) / calib.npmt;
        }
        double sqq = 0.0;
        double sqt = 0.0;
        double sqhit = 0.0;
        for (PmtTable::const_iterator it = m_pmtTable.begin(); it != m_pmtTable.end(); ++it) {
            if (!it->used) continue;
            if ( (it->type & PmtType::PMT_20INCH) != it->type ) continue;
            sqq += (it->q - calib.meanq) * (it->q - calib.meanq);
            sqt += (it->fht - calib.meant) * (it->fht - calib.meant);
            sqhit += (static_cast<double>(it->hittime.size()) - calib.meanhit) * (static_cast<double>(it->hittime.size()) - calib.meanhit);
        }
        if (calib.npmt > 1) {
            calib.stdq = std::sqrt(sqq / (calib.npmt - 1));
            calib.stdt = std::sqrt(sqt / (calib.npmt - 1));
            calib.stdhit = std::sqrt(sqhit / (calib.npmt - 1));
        }

        LogInfo << "TotQ: CD = " << calib.totq << ", WP = " << totq_wp << '\n';

        bool is_possibly_cd_muon = false;
        bool is_possibly_wp_muon = false;

        if (
            calib.totq >= m_cd_muon_totq_thold && 
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
            calib.totq < m_cd_muon_totq_thold && 
            totq_wp >= m_wp_muon_totq_thold &&
            curts - m_wp_last_muon > m_wp_afterpulse_thold
        ) {
            m_wp_last_muon = curts;
            is_possibly_wp_muon = true;
        }
        else if (
            calib.totq >= m_cd_muon_totq_thold && 
            totq_wp < m_wp_muon_totq_thold &&
            (curts - m_cd_last_muon > TimeStamp{0, 2000000} || curts - m_wp_last_muon > TimeStamp{0, 2000000})
        ) {
            m_cd_last_muon = curts;
            is_possibly_cd_muon = true;
        }

        LogInfo << "Is possibly CD muon: " << is_possibly_cd_muon << ", is possibly WP muon: " << is_possibly_wp_muon << '\n';

        std::vector<track> tracks;
        m_trkSaver.reset();
        if (is_possibly_cd_muon || is_possibly_wp_muon) {
            for (JM::NavBuffer::Iterator it = bufwrap.begin(); it != bufwrap.end(); ++it) {
                TimeStamp otherts = it->get()->TimeStamp().GetTimeSpec();
                if (curts - otherts < TimeStamp{0, -500} || TimeStamp{0, 500} < curts - otherts) continue;
                LogInfo << "Current detector: " << bufwrap.curEvt()->getDetectorType() << ", other detector: " << it->get()->getDetectorType() << '\n';
                if (is_possibly_cd_muon && it->get()->getDetectorType() == JM::EvtNavigator::DetectorType::CD) {
                    JM::CdTrackRecHeader* basic_cdt_hdr = JM::getHeaderObject<JM::CdTrackRecHeader>(bufwrap.curEvt());
                    addTrack(basic_cdt_hdr, "CdBasic", curts, tracks);
                    LogInfo << "CdBasic: " << basic_cdt_hdr << '\n';
                    JM::CdTrackRecHeader* classify_cdt_hdr = JM::getHeaderObject<JM::CdTrackRecHeader>(bufwrap.curEvt(), "/Event/CdTrackRecClassify");
                    addTrack(classify_cdt_hdr, "CdClassify", curts, tracks);
                    LogInfo << "CdClassify: " << classify_cdt_hdr << '\n';
                    RecTrks rtrks;
                    if (!m_recTool->reconstruct(&rtrks)) {
                        LogWarn << "Could not reconstruct the event with reconstruction tool\n";
                    }
                    addTrack(rtrks, "CdWpTtChi2", curts, track::loc::cd, tracks);
                    LogInfo << "CdWpTtChi2: " << rtrks.size() << '\n';
                    m_trkSaver.add(rtrks, "CdWpTtChi2", bufwrap.curEvt()->RunID(), curts);
                    m_trkSaver.add(classify_cdt_hdr, "CdClassify", bufwrap.curEvt()->RunID(), curts);
                }
                if (is_possibly_wp_muon && it->get()->getDetectorType() == JM::EvtNavigator::DetectorType::WP) {
                    JM::WpRecHeader* basic_wpt_hdr = JM::getHeaderObject<JM::WpRecHeader>(bufwrap.curEvt());
                    addTrack(basic_wpt_hdr, "WpBasic", curts, tracks);
                    LogInfo << "WpBasic: " << basic_wpt_hdr << '\n';
                    JM::WpRecHeader* classify_wpt_hdr = JM::getHeaderObject<JM::WpRecHeader>(bufwrap.curEvt(), "/Event/WpTrackRecClassify");
                    addTrack(classify_wpt_hdr, "WpClassify", curts, tracks);
                    // TODO NOT FOR NOW: Add track saver for WpClassify
                }
            }
            if (tracks.empty()) {
                tracks.push_back(track{"CdBasic", vec3{0.0, 0.0, 20000.0}, vec3{0.0, 0.0, -20000.0}, 0.0, curts, track::loc::cd, -1.0});
            }
        }
        m_trkSaver.fill();

        addTtToTrack(tracks, curts);

        if (m_tsEvt <= curts) {
            addFeature(tracks, curts, runId);
        }

        std::vector<vertex> vertices;
        // JM::OecHeader* oec_hdr = JM::getHeaderObject<JM::OecHeader>(bufwrap.curEvt());
        // addVertex(oec_hdr, "Oec", curts, totq_cd, vertices);
        // JM::CdVertexRecHeader* basic_cdv_hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(bufwrap.curEvt());
        // addVertex(basic_cdv_hdr, "Basic", curts, totq_cd, vertices);
        // JM::CdVertexRecHeader* jvertex_cdv_hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(bufwrap.curEvt(), "/Event/CdVertexRecJVertex");
        // addVertex(jvertex_cdv_hdr, "JVertex", curts, totq_cd, vertices);
        JM::CdVertexRecHeader* mixedphase_cdv_hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(bufwrap.curEvt(), "/Event/CdVertexRecMixedPhase");
        addVertex(mixedphase_cdv_hdr, "MixedPhase", curts, calib, vertices);
        JM::CdVertexRecHeader* omilrec_cdv_hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(bufwrap.curEvt(), "/Event/CdVertexRecOMILREC");
        addVertex(omilrec_cdv_hdr, "OMILREC", curts, calib, vertices);
        JM::CdVertexRecHeader* omilrec_jvertex_cdv_hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(bufwrap.curEvt(), "/Event/CdVertexRecOMILREC_JVtx");
        addVertex(omilrec_jvertex_cdv_hdr, "OMILREC_JVtx", curts, calib, vertices);

        evt->tracks = tracks;
        evt->vertices = vertices;
        EventCache::insert(curts, evt);
        LogInfo << *evt << '\n';
    }

    // DEBUG --- Timing
    auto t_after_load = clock::now();
    // DEBUG --- Timing

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

bool AnalysisGroupC::finalize() {
    if (m_loader && !m_loader->finalize()) return false;
    if (m_recTool && !(dynamic_cast<ToolBase*>(m_recTool))->finalize()) return false;
    
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