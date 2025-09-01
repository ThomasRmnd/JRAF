#include "loader/cd/CdLRangeFiller.hpp"

#include "Event/CdLpmtCalibHeader.h"
#include "EvtNavigator/EvtNavHelper.h"

#include "utils/LoadHelper.hpp"

CdLRangeFiller::CdLRangeFiller(const std::string& name, double res_20inch) :
    RangeFiller<CdGeom>{name}, 
    m_res_20inch{res_20inch}
{}

std::size_t CdLRangeFiller::size() {
    if (!m_geom) return 0ul;
    return m_geom->findPmt20inchNum();
}

bool CdLRangeFiller::loadGeom() {
    return details::loadCdGeom(m_first, m_last, m_geom, m_res_20inch, 0.0);
}

bool CdLRangeFiller::fill(JM::EvtNavigator* nav) {
    JM::CdLpmtCalibHeader* cdl_hdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
    if (cdl_hdr && cdl_hdr->event()) {
        const std::list<JM::CalibPmtChannel*>& clb_list = cdl_hdr->event()->calibPMTCol();
        LogInfo << "Loading " << clb_list.size() << " 20-inch PMTs for CD\n";
        if (!loadCol(clb_list)) return false;
    }
    return true;
}

bool CdLRangeFiller::loadCol(const std::list<JM::CalibPmtChannel*>& clb_list) {
    return details::loadCol(clb_list.begin(), clb_list.end(), m_first, m_last);
}