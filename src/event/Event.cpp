#include "event/Event.hpp"

#include "SniperKernel/SniperLog.h"

#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdTrackRecHeader.h"
#include "Event/CdVertexRecHeader.h"
#include "Event/OecHeader.h"
#include "Event/SimHeader.h"
#include "Event/TtRecHeader.h"
#include "Event/WpRecHeader.h"
#include "EvtNavigator/EvtNavHelper.h"

bool Event::load(JM::EvtNavigator* nav) {
    tracks.clear();
    vertices.clear();
    if (!nav) {
        LogError << "Event navigator is nullptr\n";
        return false;
    }
    run_id = nav->RunID();
    ts = TimeStamp{nav->TimeStamp().GetTimeSpec()};

    JM::CdLpmtCalibHeader* calib_hdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
    calibration_context calib;
    if (calib_hdr && calib_hdr->event()) {
        const std::list<JM::CalibPmtChannel*>& chlist = calib_hdr->event()->calibPMTCol();
        for (const JM::CalibPmtChannel* ch : chlist) {
            double q = static_cast<double>(ch->nPE());
            calib.totq += q;
            for (float t : ch->time()) {
                calib.meant += static_cast<double>(t);
            }
            ++calib.npmt;
            calib.nhit += ch->time().size();
            if (q < calib.minq) calib.minq = q;
            if (q > calib.maxq) calib.maxq = q;
        }
        if (calib.npmt > 0) {
            calib.meanq = calib.totq / static_cast<double>(calib.npmt);
            calib.meant = calib.meant / static_cast<double>(calib.nhit);
            calib.meanhit = static_cast<double>(calib.nhit) / static_cast<double>(calib.npmt);
        }
        double sqq = 0.0;
        double sqt = 0.0;
        double sqhit = 0.0;
        for (const JM::CalibPmtChannel* ch : chlist) {
            double q = static_cast<double>(ch->nPE());
            sqq += (q - calib.meanq) * (q - calib.meanq);
            for (float t : ch->time()) {
                sqt += (static_cast<double>(t) - calib.meant) * (static_cast<double>(t) - calib.meant);
            }
            sqhit += (static_cast<double>(ch->time().size()) - calib.meanhit) * (static_cast<double>(ch->time().size()) - calib.meanhit);
        }
        if (calib.npmt > 1) {
            calib.stdq = std::sqrt(sqq / static_cast<double>(calib.npmt - 1ul));
            calib.stdt = std::sqrt(sqt / static_cast<double>(calib.npmt - 1ul));
            calib.stdhit = std::sqrt(sqhit / static_cast<double>(calib.npmt - 1ul));
        }
    }
    
    JM::SimHeader* sim_hdr = JM::getHeaderObject<JM::SimHeader>(nav);
    if (sim_hdr) type = sim_hdr->getEventType();
    else type = "Unknown";

    // loadCdTrack(nav);
    // loadWpTrack(nav);
    // loadTtTrack(nav);
    loadCdVertex(nav, calib);
    return true;
}

void Event::loadCdTrack(JM::EvtNavigator* nav) {
    JM::CdTrackRecHeader* hdr = JM::getHeaderObject<JM::CdTrackRecHeader>(nav);
    if (!hdr || !hdr->event()) return;
    const std::vector<JM::RecTrack*>& trks = hdr->event()->tracks();
    for (const JM::RecTrack* trk : trks) {
        loadTrack("CdBasic", trk, track::loc::cd);
    }
}

void Event::loadWpTrack(JM::EvtNavigator* nav) {
    JM::WpRecHeader* hdr = JM::getHeaderObject<JM::WpRecHeader>(nav);
    if (!hdr || !hdr->event()) return;
    const std::vector<JM::RecTrack*>& trks = hdr->event()->tracks();
    for (const JM::RecTrack* trk : trks) {
        loadTrack("WpBasic", trk, track::loc::wp);
    }
}

void Event::loadTtTrack(JM::EvtNavigator* nav) {
    JM::TtRecHeader* hdr = JM::getHeaderObject<JM::TtRecHeader>(nav);
    if (!hdr || !hdr->event()) return;
    int ntracks = hdr->event()->nTracks();
    const std::vector<float>& coeff0 = hdr->event()->Coeff0();
    const std::vector<float>& coeff1 = hdr->event()->Coeff1();
    const std::vector<float>& coeff2 = hdr->event()->Coeff2();
    const std::vector<float>& coeff3 = hdr->event()->Coeff3();
    const std::vector<float>& coeff4 = hdr->event()->Coeff4();
    const std::vector<float>& coeff5 = hdr->event()->Coeff5();
    const std::vector<float>& chi2 = hdr->event()->Chi2();
    std::vector<JM::RecTrack*> trks;
    trks.reserve(ntracks);
    for (int i = 0; i < ntracks; ++i) {
        JM::RecTrack trk;
        TVector3 ipos(coeff0[i], coeff1[i], coeff2[i]);
        TVector3 dir(coeff3[i], coeff4[i], coeff5[i]);
        TVector3 fpos = ipos - 2.0 * (ipos * dir) * dir;
        trk.setStart(CLHEP::HepLorentzVector(ipos.X(), ipos.Y(), ipos.Z()));
        trk.setEnd(CLHEP::HepLorentzVector(fpos.X(), fpos.Y(), fpos.Z()));
        trk.setPESum(0.0f);
        trk.setQuality(chi2[i]);
        loadTrack("TtBasic", &trk, track::loc::tt);
    }
}

void Event::loadCdVertex(JM::EvtNavigator* nav, const calibration_context& calib) {
    // JM::CdVertexRecHeader* hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(nav);
    // if (!hdr || !hdr->event()) return;
    // const std::vector<JM::RecVertex*>& vtxs = hdr->event()->vertices();
    // for (const JM::RecVertex* vtx : vtxs) {
    //     loadVertex(vtx);
    // }
    JM::OecHeader* hdr = JM::getHeaderObject<JM::OecHeader>(nav);
    if (!hdr || !hdr->event("JM::OecEvt")) return;
    JM::OecEvt* evt = dynamic_cast<JM::OecEvt*>(hdr->event("JM::OecEvt"));
    vertices.emplace_back(
        "Oec",
        vec3{evt->getVertexX(), evt->getVertexY(), evt->getVertexZ()},
        evt->getEnergy(), ts, calib, type
    );
}

void Event::loadTrack(const std::string& method, const JM::RecTrack* trk, const track::loc& det) {
    if (!trk) return;
    tracks.emplace_back(method, *trk, ts, det);
}