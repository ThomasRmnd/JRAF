#include "event/EventBuilder.hpp"

#include "SniperKernel/ToolFactory.h"

#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdSpmtCalibHeader.h"
#include "Event/TtCalibHeader.h"
#include "Event/WpCalibHeader.h"
#include "EvtNavigator/EvtNavHelper.h"
#include "UtilsThomas/utils/DetectorType.hpp"

#include "event/EventCache.hpp"
#include "utils/NavBufferWrapper.hpp"

DECLARE_TOOL(EventBuilder);

EventBuilder::EventBuilder(const std::string& name) : 
    ToolBase{name} 
{
    declProp("MuonTagger", m_mutagger_name = "MuonTagger");

    declProp("Loader", m_loader_name = "JointLoader");
    declProp("CdFiller", m_cdfiller_name = "CdRangeFiller");
    declProp("WpFiller", m_wpfiller_name = "WpRangeFiller");
    declProp("TtFiller", m_ttfiller_name = "TtRangeFiller");

    declProp("RecTool", m_rectool_name);
    declProp("TtReco", m_tt_reco_file.filename = "");
}

bool EventBuilder::initialize() {
    AlgBase* alg = m_par->findAlg("JRAF");
    if (!alg) {
        LogError << "Cannot retrieve parent algorithm\n";
        return false;
    }

    m_mutagger = alg->tool<MuonTagger>(m_mutagger_name);
    if (!m_mutagger) {
        LogError << "Failed to retrieve MuonTagger tool named " << m_mutagger_name << '\n';
        return false;
    }
    if (!m_mutagger->initialize()) return false;


    m_loader = alg->tool<Loader>(m_loader_name);
    if (!m_loader) {
        LogError << "Failed to retrieve reconstruction tool named " << m_loader_name << '\n';
        return false;
    }
    RangeFiller<CdFillerTag>* cd_filler = alg->tool<RangeFiller<CdFillerTag>>(m_cdfiller_name);
    RangeFiller<WpFillerTag>* wp_filler = alg->tool<RangeFiller<WpFillerTag>>(m_wpfiller_name);
    RangeFiller<TtFillerTag>* tt_filler = alg->tool<RangeFiller<TtFillerTag>>(m_ttfiller_name);
    if (!m_loader->configure(&m_pmt_table, cd_filler, wp_filler, tt_filler)) return false;
	if (!m_loader->initialize()) return false;

    m_rectool = alg->tool<IRecMuonTool>(m_rectool_name);
    if (!m_rectool) {
        LogError << "Failed to retrieve reconstruction tool named " << m_rectool_name << '\n';
        return false;
    }
    if (!m_rectool->configure(&m_params, &m_pmt_table)) return false;
	// if (!dynamic_cast<ToolBase*>(m_recTool)->initialize()) return false;

    if (!m_tt_reco_file.init()) return false;

    return true;
}

bool EventBuilder::finalize() {
    if (m_mutagger && !m_mutagger->finalize()) return false;
    if (m_loader && !m_loader->finalize()) return false;
    if (m_rectool && !(dynamic_cast<ToolBase*>(m_rectool))->finalize()) return false;
    return true;
}

DetectorType EventBuilder::getDetectorType(JM::EvtNavigator* nav) {
    DetectorType type = DetectorType::UNKNOWN;

    JM::EvtNavigator::DetectorType evt_type = nav->getDetectorType();

    if (evt_type == JM::EvtNavigator::DetectorType::CD) type |= DetectorType::CD;
    if (evt_type == JM::EvtNavigator::DetectorType::WP) type |= DetectorType::WP;
    if (evt_type == JM::EvtNavigator::DetectorType::TT) type |= DetectorType::TT;

    return type;
}

calibration_context EventBuilder::getCalibrationContext(const std::list<JM::CalibPmtChannel*>& clb_list) {
    calibration_context calib;
    for (JM::CalibPmtChannel* clb : clb_list) {
        if (!clb) continue;
        for (float t : clb->time()) {
            calib.meant += static_cast<double>(t);
        }
        ++calib.npmt;
        calib.nhit += clb->time().size();
        double totq = static_cast<double>(clb->sumCharge());
        calib.totq += totq;
        if (totq < calib.minq) calib.minq = totq;
        if (totq > calib.maxq) calib.maxq = totq;
    }
    if (calib.npmt > 0) {
        calib.meanq = calib.totq / static_cast<double>(calib.npmt);
        calib.meant = calib.meant / static_cast<double>(calib.nhit);
        calib.meanhit = static_cast<double>(calib.nhit) / static_cast<double>(calib.npmt);
    }
    double sqq = 0.0;
    double sqt = 0.0;
    double sqhit = 0.0;
    for (JM::CalibPmtChannel* clb : clb_list) {
        if (!clb) continue;
        double totq = static_cast<double>(clb->sumCharge());
        sqq += (totq - calib.meanq) * (totq - calib.meanq);
        for (float t : clb->time()) {
            sqt += (static_cast<double>(t) - calib.meant) * (static_cast<double>(t) - calib.meant);
        }
        sqhit += (static_cast<double>(clb->time().size()) - calib.meanhit) * (static_cast<double>(clb->time().size()) - calib.meanhit);
    }
    if (calib.npmt > 1) {
        calib.stdq = std::sqrt(sqq / static_cast<double>(calib.npmt - 1ul));
        calib.stdt = std::sqrt(sqt / static_cast<double>(calib.nhit - 1ul));
        calib.stdhit = std::sqrt(sqhit / static_cast<double>(calib.npmt - 1ul));
    }
    return calib;
}

void EventBuilder::addTrack(RecTrks& rec_tracks, const std::string& method, double totq_cd, double totq_wp, const TimeStamp& ts, const track::loc& det, std::vector<track>& tracks) {
    for (int k = 0; k < rec_tracks.size(); ++k) {
        tracks.push_back(track{
            method, vec3{rec_tracks.getStart(k)}, vec3{rec_tracks.getEnd(k)}, totq_cd, totq_wp, ts, det, rec_tracks.getQuality(k)
        });
    }
}

void EventBuilder::addTrack(JM::CdTrackRecHeader* cdt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks) {
    if (!cdt_hdr || !cdt_hdr->event()) return;
    const std::vector<JM::RecTrack*>& rec_tracks = cdt_hdr->event()->tracks();
    for (JM::RecTrack* t : rec_tracks) {
        tracks.push_back(track{
            method, *t, ts, track::loc::cd
        });
    }
}

void EventBuilder::addTrack(JM::WpRecHeader* wpt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks) {
    if (!wpt_hdr || !wpt_hdr->event()) return;
    const std::vector<JM::RecTrack*>& rec_tracks = wpt_hdr->event()->tracks();
    for (JM::RecTrack* t : rec_tracks) {
        tracks.push_back(track{
            method, *t, ts, track::loc::wp
        });
    }
}

void EventBuilder::addTrack(JM::TtRecHeader* ttt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks) {
    if (!ttt_hdr || !ttt_hdr->event()) return;
    JM::TtRecEvt* ttt_evt = ttt_hdr->event();
    for (int k = 0; k < ttt_evt->nTracks(); ++k) {
        vec3 ipos{ttt_evt->Coeff0()[k], ttt_evt->Coeff1()[k], ttt_evt->Coeff2()[k]};
        vec3 dir = unit(vec3{ttt_evt->Coeff3()[k], ttt_evt->Coeff4()[k], ttt_evt->Coeff5()[k]});
        vec3 fpos = ipos - 2.0 * dot(ipos, dir) * dir;
        tracks.push_back(track{
            method, ipos, fpos, 0.0, 0.0, ts, track::loc::tt, ttt_evt->Chi2()[k]
        });
    }
}

void EventBuilder::addVertex(JM::OecHeader* oec_hdr, const std::string& method, const TimeStamp& ts, const calibration_context& calib, std::vector<vertex>& vertices) {
    if (!oec_hdr || !oec_hdr->event("JM::OecEvt")) return;
    JM::OecEvt* oec_evt = dynamic_cast<JM::OecEvt*>(oec_hdr->event("JM::OecEvt"));
    vertices.push_back(vertex{
        method, vec3{oec_evt->getVertexX(), oec_evt->getVertexY(), oec_evt->getVertexZ()}, oec_evt->getEnergy(), ts, calib, "Unknown"
    });
}

void EventBuilder::addVertex(JM::CdVertexRecHeader* cdv_hdr, const std::string& method, const TimeStamp& ts, const calibration_context& calib, std::vector<vertex>& vertices) {
    if (!cdv_hdr || !cdv_hdr->event()) return;
    const std::vector<JM::RecVertex*>& rec_vertices = cdv_hdr->event()->vertices();
    for (JM::RecVertex* v : rec_vertices) {
        vertices.push_back(vertex{
            method, vec3{v->x(), v->y(), v->z()}, v->energy(), ts, calib, "Unknown"
        });
    }
}

int EventBuilder::getTtLayerId(double z) {
    if (24000.0 <= z && z <= 25000.0) return 0;  // main
    if (25500.0 <= z && z <= 26500.0) return 1;  // main
    if (27000.0 <= z && z <= 28000.0) return 2;  // main
    if (30000.0 <= z && z <= 30200.0) return 3;  // chimney
    if (30200.0 <= z && z <= 30400.0) return 4;  // chimney
    if (30400.0 <= z && z <= 30600.0) return 5;  // chimney

    return -1; // not inside any valid layer
};

void EventBuilder::addTtToTrack(std::vector<track>& tracks, const TimeStamp& curts) {
    if (!m_tt_reco_file.find(curts)) return;
    
    if (m_tt_reco_file.NTracks != 1) {
        LogInfo << "No info or bundle muons reconstructed by the TT\n";
        return;
    }
    if (m_tt_reco_file.NPoints[0] < 3) {
        LogInfo << "Track has less than 3 points\n";
        return;
    }
    std::unordered_set<int> layers_hit;
    layers_hit.reserve(6);
    for (int i = 0; i < m_tt_reco_file.NTotPoints; ++i) {
        int lid = getTtLayerId(m_tt_reco_file.PointZ[i] + 26452.0);
        if (lid < 0) continue;
        layers_hit.insert(lid);
    }
    if (layers_hit.size() < 3) {
        LogInfo << "Track is not in three different layers of the TT\n";
        return;
    }

    vec3 ipos{m_tt_reco_file.Coeff0[0], m_tt_reco_file.Coeff1[0], m_tt_reco_file.Coeff2[0] + 26452.0};
    vec3 dir = unit(vec3{m_tt_reco_file.Coeff3[0], m_tt_reco_file.Coeff4[0], m_tt_reco_file.Coeff5[0]});
    vec3 fpos = ipos - 2.0 * dot(ipos, dir) * dir;
    tracks.push_back(track{
        "Tt", ipos, fpos, 0.0, 0.0, curts, track::loc::tt, m_tt_reco_file.Chi2[0]
    });
}

bool EventBuilder::build(JM::NavBuffer* buf) {
    for (NavBufferWrapper bufwrap(*buf); bufwrap.current() != bufwrap.end(); bufwrap.next()) {
        if (EventCache::contains(bufwrap.curEvt())) continue;

        JM::EvtNavigator* curnav = bufwrap.curEvt();
        if (!curnav) {
            LogError << "EvtNavigator is nullptr\n";
            return false;
        }

        std::shared_ptr<Event> evt = std::make_shared<Event>();

        TimeStamp curts{curnav->TimeStamp().GetTimeSpec()};
        DetectorType curdet = getDetectorType(curnav);
        if (curdet == DetectorType::UNKNOWN) {
            LogError << "Unknown detector type\n";
            return false;
        }

        calibration_context calib_cd, calib_wp;

        JM::CdLpmtCalibHeader* cdl_calib_hdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(curnav);
        JM::EvtNavigator* cdl_evt_nav = nullptr;
        if (cdl_calib_hdr && cdl_calib_hdr->event()) {
            calib_cd = getCalibrationContext(cdl_calib_hdr->event()->calibPMTCol());
            cdl_evt_nav = curnav;
        }
        JM::WpCalibHeader* wp_calib_hdr = JM::getHeaderObject<JM::WpCalibHeader>(curnav);
        JM::EvtNavigator* wp_evt_nav = nullptr;
        if (wp_calib_hdr && wp_calib_hdr->event()) {
            calib_wp = getCalibrationContext(wp_calib_hdr->event()->calibPMTCol());
            wp_evt_nav = curnav;
        }

        DetectorType jointdet = curdet;
        for (JM::NavBuffer::Iterator it = bufwrap.begin(); it != bufwrap.end(); ++it) {
            if (it == bufwrap.current()) continue;
            JM::EvtNavigator* othernav = it->get();
            if (!othernav) {
                LogError << "EvtNavigator is nullptr\n";
                return false;
            }
            TimeStamp otherts{othernav->TimeStamp().GetTimeSpec()};
            if (curts - otherts < TimeStamp{0, -500} || TimeStamp{0, 500} < curts - otherts) continue;
            DetectorType otherdet = getDetectorType(othernav);
            if (otherdet == DetectorType::UNKNOWN) {
                LogError << "Unknown detector type\n";
                return false;
            }
            if ( (jointdet & otherdet) != DetectorType::UNKNOWN ) continue;
            cdl_calib_hdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(othernav);
            if (cdl_calib_hdr && cdl_calib_hdr->event()) {
                calib_cd = getCalibrationContext(cdl_calib_hdr->event()->calibPMTCol());
                cdl_evt_nav = othernav;
            }
            wp_calib_hdr = JM::getHeaderObject<JM::WpCalibHeader>(othernav);
            if (wp_calib_hdr && wp_calib_hdr->event()) {
                calib_wp = getCalibrationContext(wp_calib_hdr->event()->calibPMTCol());
                wp_evt_nav = othernav;
            }
            // Could change the reference time for PMT hits but not necessary here (only necessary for joint loader)
            jointdet |= otherdet;
        }

        LogInfo << "TotQ: CD = " << calib_cd.totq << ", WP = " << calib_wp.totq << '\n';
        
        evt->run_id = curnav->RunID();
        evt->ts = curts;

        muon_tagging_context muctx{
            curts, curdet,
            calib_cd, calib_wp,
            cdl_evt_nav, wp_evt_nav
        };
        muon_tagging_result mures = m_mutagger->tag(muctx);

        if (mures.is_possibly_cd_muon && mures.is_possibly_wp_muon) {
            evt->totq_cd = calib_cd.totq;
            evt->totq_wp = calib_wp.totq;
        }
        else if (mures.is_possibly_wp_muon) {
            evt->totq_cd = 0.0;
            evt->totq_wp = calib_wp.totq;
        }
        else if (mures.is_possibly_cd_muon) {
            evt->totq_cd = calib_cd.totq;
            evt->totq_wp = 0.0;
        }

        LogInfo << "Is possibly CD muon: " << mures.is_possibly_cd_muon << ", is possibly WP muon: " << mures.is_possibly_wp_muon << '\n';

        std::vector<track> tracks;
        if (mures.is_possibly_cd_muon || mures.is_possibly_wp_muon) {
            LoadingResult loadres = m_loader->load(&bufwrap);
            if (!loadres.ok) return false;

            if (cdl_evt_nav) {
                JM::CdTrackRecHeader* basic_cdt_hdr = JM::getHeaderObject<JM::CdTrackRecHeader>(cdl_evt_nav);
                addTrack(basic_cdt_hdr, "CdBasic", curts, tracks);
                LogInfo << "CdBasic: " << basic_cdt_hdr << '\n';
                JM::CdTrackRecHeader* classify_cdt_hdr = JM::getHeaderObject<JM::CdTrackRecHeader>(cdl_evt_nav, "/Event/CdTrackRecClassify");
                addTrack(classify_cdt_hdr, "CdClassify", curts, tracks);
                LogInfo << "CdClassify: " << classify_cdt_hdr << '\n';
                RecTrks rtrks;
                if (!m_rectool->reconstruct(&rtrks)) {
                    LogWarn << "Could not reconstruct the event with reconstruction tool\n";
                }
                addTrack(rtrks, "CdWpTtChi2", calib_cd.totq, calib_wp.totq, curts, track::loc::cd, tracks);
                LogInfo << "CdWpTtChi2: " << rtrks.size() << '\n';
            }
            if (wp_evt_nav) {
                JM::WpRecHeader* basic_wpt_hdr = JM::getHeaderObject<JM::WpRecHeader>(wp_evt_nav);
                LogInfo << "WpBasic: " << basic_wpt_hdr << '\n';
                addTrack(basic_wpt_hdr, "WpBasic", curts, tracks);
                JM::WpRecHeader* classify_wpt_hdr = JM::getHeaderObject<JM::WpRecHeader>(wp_evt_nav, "/Event/WpTrackRecClassify");
                LogInfo << "WpClassify: " << classify_wpt_hdr << '\n';
                addTrack(classify_wpt_hdr, "WpClassify", curts, tracks);
                // TODO NOT FOR NOW: Add track saver for WpClassify
            }
            
            addTtToTrack(tracks, curts);
            
            if (tracks.empty()) {
                tracks.push_back(track{"Default", vec3{0.0, 0.0, 20000.0}, vec3{0.0, 0.0, -20000.0}, calib_cd.totq, calib_wp.totq, curts, track::loc::cd, -1.0});
            }
        }

        std::vector<vertex> vertices;
        JM::CdVertexRecHeader* mixedphase_cdv_hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(curnav, "/Event/CdVertexRecMixedPhase");
        addVertex(mixedphase_cdv_hdr, "MixedPhase", curts, calib_cd, vertices);
        JM::CdVertexRecHeader* omilrec_cdv_hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(curnav, "/Event/CdVertexRecOMILREC");
        addVertex(omilrec_cdv_hdr, "OMILREC", curts, calib_cd, vertices);
        JM::CdVertexRecHeader* omilrec_jvertex_cdv_hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(curnav, "/Event/CdVertexRecOMILREC_JVtx");
        addVertex(omilrec_jvertex_cdv_hdr, "OMILREC_JVtx", curts, calib_cd, vertices);

        evt->tracks = tracks;
        evt->vertices = vertices;
        EventCache::insert(curts, evt);
        LogInfo << *evt << '\n';
    }

    return true;
}