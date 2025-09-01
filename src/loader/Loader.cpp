#include "loader/Loader.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperPtr.h"

#include "Geometry/IRecGeomSvc.hh"

Loader::Loader(
    const std::string& name, 
    PmtTable* pmt_table, 
    const std::shared_ptr<RangeFiller<CdGeom>>& cd_filler, 
    const std::shared_ptr<RangeFiller<WpGeom>>& wp_filler, 
    const std::shared_ptr<RangeFiller<TtGeom>>& tt_filler, IRecGeomSvc* geom_svc
) :
    ToolBase{name},
    m_pmt_table{pmt_table},
    m_geom_svc{geom_svc},
    m_cd_filler{cd_filler},
    m_wp_filler{wp_filler},
    m_tt_filler{tt_filler}
{}

bool Loader::initialize() {
    if (!initGeomSvc()) return false;
    if (m_cd_filler) {
        m_cd_filler->setRange(m_cd_it, m_wp_it);
        if (!m_cd_filler->initialize()) return false;
    }
    if (m_wp_filler) {
        m_wp_filler->setRange(m_wp_it, m_tt_it);
        if (!m_wp_filler->initialize()) return false;
    }
    if (m_tt_filler) {
        m_tt_filler->setRange(m_tt_it, m_pmt_table->end());
        if (!m_tt_filler->initialize()) return false;
    }
    return true;
}

bool Loader::initGeomSvc() {
    if (!m_geom_svc) {
        LogError << "IRecGeomSvc has not been set!\n";
        return false;
    }
    
    std::size_t n_cd = 0ul, n_wp = 0ul, n_tt = 0ul;
    if (m_cd_filler) {
        m_cd_filler->setGeom(m_geom_svc->getCdGeom());
        n_cd = m_cd_filler->size();
    }
    if (m_wp_filler) {
        m_wp_filler->setGeom(m_geom_svc->getWpGeom());
        n_wp = m_wp_filler->size();
    }
    if (m_tt_filler) {
        m_tt_filler->setGeom(m_geom_svc->getTtGeom());
        n_tt = m_tt_filler->size();
    }

    LogInfo << "Currently loading: " << n_cd << " CD PMTs, " << n_wp << " WP PMTs, " << n_tt << " TT channels\n";

    m_pmt_table->resize(n_cd + n_wp + n_tt);

    m_cd_it = m_pmt_table->begin(); // always at beginning whether CD is loaded or not
    m_wp_it = m_cd_it + n_cd;
    m_tt_it = m_wp_it + n_wp;

    LogDebug << "table.begin = " << &(*(m_pmt_table->begin())) << 
                ", cd_it = " << &(*m_cd_it) << 
                ", wp_it = " << &(*m_wp_it) << 
                ", tt_it = " << &(*m_tt_it) << 
                ", table.end = " << &(*(m_pmt_table->end())) << std::endl;

    return true;
}

bool Loader::finalize() {
    if (m_cd_filler && !m_cd_filler->finalize()) return false;
    if (m_wp_filler && !m_wp_filler->finalize()) return false;
    if (m_tt_filler && !m_tt_filler->finalize()) return false;
    return true;
}

bool Loader::loadNav(JM::EvtNavigator* nav) {
    if (m_cd_filler && !m_cd_filler->fill(nav)) return false;
    if (m_wp_filler && !m_wp_filler->fill(nav)) return false;
    if (m_tt_filler && !m_tt_filler->fill(nav)) return false;
    return true;
}

void Loader::unload() {
    for (PmtProp& pmt : *m_pmt_table) {
        pmt.used = false;
    }
}