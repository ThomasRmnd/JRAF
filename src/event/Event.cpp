#include "event/Event.hpp"

#include "SniperKernel/SniperLog.h"

#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdTrackRecHeader.h"
#include "Event/CdVertexRecHeader.h"
#include "Event/SimHeader.h"
#include "Event/WpRecHeader.h"
#include "EvtNavigator/EvtNavHelper.h"

bool Event::load(JM::EvtNavigator* nav) {
    tracks.clear();
    vertices.clear();
    if (!nav) {
        LogError << "Event navigator is nullptr\n";
        return false;
    }
    ts = TimeStamp{nav->TimeStamp().GetTimeSpec()};

    JM::CdLpmtCalibHeader* calib_hdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
    totq = 0.0;
    if (calib_hdr && calib_hdr->event()) {
        const std::list<JM::CalibPmtChannel*>& chlist = calib_hdr->event()->calibPMTCol();
        for (const JM::CalibPmtChannel* ch : chlist) {
            totq += ch->nPE();
        }
    }
    
    JM::SimHeader* sim_hdr = JM::getHeaderObject<JM::SimHeader>(nav);
    if (sim_hdr) type = sim_hdr->getEventType();
    else type = "Unknown";

    loadCdTrack(nav);
    loadWpTrack(nav);
    loadCdVertex(nav);
    return true;
}

void Event::loadCdTrack(JM::EvtNavigator* nav) {
    JM::CdTrackRecHeader* hdr = JM::getHeaderObject<JM::CdTrackRecHeader>(nav);
    if (!hdr || !hdr->event()) return;
    const std::vector<JM::RecTrack*>& trks = hdr->event()->tracks();
    for (const JM::RecTrack* trk : trks) {
        loadTrack(trk, track::loc::cd);
    }
}

void Event::loadWpTrack(JM::EvtNavigator* nav) {
    JM::WpRecHeader* hdr = JM::getHeaderObject<JM::WpRecHeader>(nav);
    if (!hdr || !hdr->event()) return;
    const std::vector<JM::RecTrack*>& trks = hdr->event()->tracks();
    for (const JM::RecTrack* trk : trks) {
        loadTrack(trk, track::loc::wp);
    }
}

void Event::loadCdVertex(JM::EvtNavigator* nav) {
    JM::CdVertexRecHeader* hdr = JM::getHeaderObject<JM::CdVertexRecHeader>(nav);
    if (!hdr || !hdr->event()) return;
    const std::vector<JM::RecVertex*>& vtxs = hdr->event()->vertices();
    for (const JM::RecVertex* vtx : vtxs) {
        loadVertex(vtx);
    }
}

void Event::loadTrack(const JM::RecTrack* trk, const track::loc& det) {
    if (!trk) return;
    tracks.emplace_back(*trk, ts, det);
}

void Event::loadVertex(const JM::RecVertex* vtx) {
    if (!vtx) return;
    vertices.emplace_back(*vtx, totq, ts, type);
}