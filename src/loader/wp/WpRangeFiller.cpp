#include "loader/wp/WpRangeFiller.hpp"

#include "Event/WpCalibHeader.h"
#include "EvtNavigator/EvtNavHelper.h"
#include "Identifier/WpID.h"

#include "utils/LoadHelper.hpp"

WpRangeFiller::WpRangeFiller(const std::string& name, double res_wp) : 
    RangeFiller<WpGeom>{name}, 
    m_res_wp{res_wp} 
{}

std::size_t WpRangeFiller::size() {
    if (!m_geom) return 0ul;
    return m_geom->getPmtNum();
}

bool WpRangeFiller::loadGeom() {
    for (PmtTable::iterator it = m_first; it != m_last; ++it) {
        Identifier id = Identifier(WpID::id(WpID::offset_WP_VPMT(std::distance(m_first, it)), 0));
        PmtGeom* pmt = m_geom->getPmt(id);
        if (!pmt) {
            LogError << "Failed to get PMT for WP of ID: " << id << '\n';
            return false;
        }
        it->pmtid = id.getValue();
        it->value = static_cast<int>(id.getValue());
        it->used = false;
        it->pos = pmt->getCenter();
        it->loc = 2;
        it->res = m_res_wp;
        it->type = PmtType::PMT_WP;
    }
    return true;
}

bool WpRangeFiller::fill(JM::EvtNavigator* nav) {
    JM::WpCalibHeader* wp_hdr = JM::getHeaderObject<JM::WpCalibHeader>(nav);
    if (wp_hdr && wp_hdr->event()) {
        const std::list<JM::CalibPmtChannel*>& clb_list = wp_hdr->event()->calibPMTCol();
        LogInfo << "Loading " << clb_list.size() << " WP PMTs\n";
        if (!loadCol(clb_list)) return false;
    }
    return true;
}

bool WpRangeFiller::loadCol(const std::list<JM::CalibPmtChannel*>& clb_list) {
    return details::loadCol(clb_list.begin(), clb_list.end(), m_first, m_last);
}