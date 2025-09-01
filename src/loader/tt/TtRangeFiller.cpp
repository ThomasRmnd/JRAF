#include "loader/tt/TtRangeFiller.hpp"

#include "Event/TtCalibHeader.h"
#include "EvtNavigator/EvtNavHelper.h"
#include "Identifier/TtID.h"
#include "Identifier/JunoDetectorID.h"

TtRangeFiller::TtRangeFiller(const std::string& name, double res_tt, ITTGeomSvc* tt_geom) :
    RangeFiller<TtGeom>{name}, 
    m_res_tt{res_tt}, 
    m_tt_geom{tt_geom} 
{}

bool TtRangeFiller::initialize() {
    if (!initTTGeomSvc()) return false;
    return RangeFiller<TtGeom>::initialize();
}

bool TtRangeFiller::initTTGeomSvc() {
    if (!m_tt_geom) {
        LogError << "TTGeomSvc has not been set\n";
        return false;
    }
    return true;
}

std::size_t TtRangeFiller::size() {
    if (!m_geom) return 0ul;
    return m_geom->getChannelNum();
}

bool TtRangeFiller::loadGeom() {
    std::map<Identifier, PmtGeom*>::const_iterator itgeom = m_geom->cbegin();
    for (PmtTable::iterator it = m_first; it != m_last; ++it) {
        if (itgeom == m_geom->cend()) {
            LogError << "Failed to get PmtGeom instance of index: " << std::distance(m_first, it) << '\n';
            return false;
        }
        it->pmtid = itgeom->first;
        it->value = static_cast<int>(itgeom->first);
        it->used = false;
        it->pos = itgeom->second->getCenter();
        it->loc = 3;
        it->res = m_res_tt;
        it->type = Pmttype::_PMTTT;
        ++itgeom;
    }
    return true;
}

bool TtRangeFiller::fill(JM::EvtNavigator* nav) {
    JM::TtCalibHeader* tt_hdr = JM::getHeaderObject<JM::TtCalibHeader>(nav);
    if (tt_hdr && tt_hdr->event()) {
        const std::list<JM::CalibTtChannel*>& clb_list = tt_hdr->event()->calibTtCol();
        LogInfo << "Loading " << clb_list.size() << " TT PMTs\n";
        if (!loadCol(clb_list)) return false;
    }
    return true;
}

bool TtRangeFiller::loadCol(const std::list<JM::CalibTtChannel*>& clb_list) {
    PmtTable::iterator it = m_first;
    for (const JM::CalibTtChannel* pmt : clb_list) {
        Identifier id = TtID::id(pmt->pmtId());
        if (std::distance(m_first, m_last) <= std::distance(m_first, it)) {
            LogError << "Mismatch between data and geometry, ID: " << id << ", index: " << std::distance(m_first, it) << '\n';
        }
        it->pmtid = pmt->pmtId();
        it->value = static_cast<int>(pmt->pmtId());
        it->used = true;
        it->q = pmt->nPE();
        it->fht = pmt->time();
        it->pos = TVector3(
            m_tt_geom->getChannelPos(id, 0), 
            m_tt_geom->getChannelPos(id, 1), 
            m_geom->getChannel(id)->getCenter().Z()
        );
        ++it;
    }
    return true;
}