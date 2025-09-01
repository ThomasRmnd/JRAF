#include "loader/cd/CdSRangeFiller.hpp"

#include "Event/CdSpmtCalibHeader.h"
#include "EvtNavigator/EvtNavHelper.h"
#include "Identifier/CdID.h"

CdSRangeFiller::CdSRangeFiller(const std::string& name, double res_3inch) :
    RangeFiller<CdGeom>{name},
    m_res_3inch{res_3inch}
{}

std::size_t CdSRangeFiller::size() {
    if (!m_geom) return 0ul;
    return m_geom->findPmt3inchNum();
}

bool CdSRangeFiller::loadGeom() {
    std::ptrdiff_t nb_lpmts = static_cast<std::ptrdiff_t>(m_geom->findPmt20inchNum());
    for (PmtTable::iterator it = m_first; it != m_last; ++it) {
        Identifier id = Identifier(CdID::id(CdID::offset_cd_spmt(nb_lpmts + std::distance(m_first, it)), 0));
        if (!CdID::is3inch(id)) {
            LogError << "PMT (id: " << id << ") is not a 3-inch PMT\n";
            return false;
        }
        PmtGeom* pmt = m_geom->getPmt(id);
        if (!pmt) {
            LogError << "Failed to get PMT (id: " << id << ") for CD\n";
            return false;
        }
        it->pmtid = id.getValue();
        it->value = id.getValue();
        it->used = false;
        it->pos = pmt->getCenter();
        it->loc = 1;
        it->res = m_res_3inch;
        it->type = Pmttype::_PMTINCH3;
    }
    return true;
}

bool CdSRangeFiller::fill(JM::EvtNavigator* nav) {
    JM::CdSpmtCalibHeader* cds_hdr = JM::getHeaderObject<JM::CdSpmtCalibHeader>(nav);
    if (cds_hdr && cds_hdr->event()) {
        const std::list<JM::CalibPmtChannel*>& clb_list = cds_hdr->event()->calibPMTCol();
        LogInfo << "Loading " << clb_list.size() << " 3-inch PMTs for CD\n";
        if (!loadCol(clb_list)) return false;
    }
    return true;
}

bool CdSRangeFiller::loadCol(const std::list<JM::CalibPmtChannel*>& clb_list) {
    int nb_lpmts = static_cast<std::ptrdiff_t>(m_geom->findPmt20inchNum());
    for (std::list<JM::CalibPmtChannel*>::const_iterator it = clb_list.begin(); it != clb_list.end(); ++it) {
        Identifier id = Identifier((*it)->pmtId());
        int i = CdID::module(id) - nb_lpmts;

        if (i < 0 || std::distance(m_first, m_last) <= i) {
            LogError << "PMT (id: " << id << ") is out of range. Index is " << i << ", while the range is " << std::distance(m_first, m_last) << '\n';
            return false;
        }

        PmtTable::iterator it_table = m_first + static_cast<std::ptrdiff_t>(i);
        it_table->used = true;
        it_table->q = (*it)->nPE();
        it_table->fht = (*it)->firstHitTime();
        it_table->hitq = (*it)->charge();
        it_table->hittime = (*it)->time();
    }
    return true;
}