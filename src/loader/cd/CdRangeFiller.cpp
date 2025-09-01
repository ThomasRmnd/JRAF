#include "loader/cd/CdRangeFiller.hpp"

#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdSpmtCalibHeader.h"
#include "EvtNavigator/EvtNavHelper.h"

#include "utils/LoadHelper.hpp"

CdRangeFiller::CdRangeFiller(const std::string& name, double res_20inch, double res_3inch) :
    RangeFiller<CdGeom>{name}, 
    m_res_20inch{res_20inch}, 
    m_res_3inch{res_3inch}
{}

std::size_t CdRangeFiller::size() {
    if (!m_geom) return 0ul;
    return m_geom->getPmtNum();
}

bool CdRangeFiller::loadGeom() {
    return details::loadCdGeom(m_first, m_last, m_geom, m_res_20inch, m_res_3inch);
}

bool CdRangeFiller::fill(JM::EvtNavigator* nav) {
    JM::CdLpmtCalibHeader* cdl_hdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
    if (cdl_hdr && cdl_hdr->event()) {
        const std::list<JM::CalibPmtChannel*>& clb_list = cdl_hdr->event()->calibPMTCol();
        LogInfo << "Loading " << clb_list.size() << " 20-inch PMTs for CD\n";
        if (!loadCol(clb_list)) return false;
    }
    JM::CdSpmtCalibHeader* cds_hdr = JM::getHeaderObject<JM::CdSpmtCalibHeader>(nav);
    if (cds_hdr && cds_hdr->event()) {
        const std::list<JM::CalibPmtChannel*>& clb_list = cds_hdr->event()->calibPMTCol();
        LogInfo << "Loading " << clb_list.size() << " 3-inch PMTs for CD\n";
        if (!loadCol(clb_list)) return false;
    }
    return true;
}

bool CdRangeFiller::loadCol(const std::list<JM::CalibPmtChannel*>& clb_list) {
    return details::loadCol(clb_list.begin(), clb_list.end(), m_first, m_last);
}