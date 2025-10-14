#include "loader/wp/AkiraWpRangeFiller.h"

#include "SniperKernel/ToolFactory.h"

#include "Event/WpCalibHeader.h"
#include "EvtNavigator/EvtNavHelper.h"
#include "Identifier/WpID.h"

#include "utils/load_helper.h"

DECLARE_TOOL(AkiraWpRangeFiller);

AkiraWpRangeFiller::AkiraWpRangeFiller(const std::string& name) : RangeFiller<WpGeom>{name} {}

AkiraWpRangeFiller::AkiraWpRangeFiller(const std::string& name, double res_wp) : 
    RangeFiller<WpGeom>{name}, 
    m_res_wp{res_wp} 
{}

std::size_t AkiraWpRangeFiller::size() {
    if (!initGeomSvc()) {
        throw std::runtime_error("You need first to set the geometry.");
    }
    return m_geom->getPmtNum();
}

bool AkiraWpRangeFiller::loadGeom() {
    for (PmtTable::iterator it = m_first; it != m_last; ++it) {
        Identifier id = Identifier(WpID::id(WpID::offset_WP_VPMT(std::distance(m_first, it)), 0));
        PmtGeom* pmt = m_geom->getPmt(id);
        if (!pmt) {
            LogError << "Failed to get PMT for WP of ID: " << id << "." << std::endl;
            return false;
        }
        it->pmtid = std::distance(m_first, it); // should be `id.getValue();` LIKE WHY IS SOMEONE USING RELATIVE COPYNO INSTEAD OF PMT ID, LIKE FY
        it->used = false;
        it->pos = pmt->getCenter();
        it->loc = 2;
        it->res = m_res_wp;
        it->type = PmtType::PMT_WP;
    }
    return true;
}

bool AkiraWpRangeFiller::fill(JM::EvtNavigator* nav) {
    JM::WpCalibHeader* wp_hdr = JM::getHeaderObject<JM::WpCalibHeader>(nav);
    if (wp_hdr) {
        if (wp_hdr->event()) {
            const std::list<JM::CalibPmtChannel*>& clb_list = wp_hdr->event()->calibPMTCol();
            LogInfo << "Loading " << clb_list.size() << " WP PMTs." << std::endl;
            if (!loadCol(clb_list)) return false;
        }
    }
    return true;
}

bool AkiraWpRangeFiller::loadCol(const std::list<JM::CalibPmtChannel*>& clb_list) {
    return details::loadCol(clb_list.begin(), clb_list.end(), m_first, m_last);
}