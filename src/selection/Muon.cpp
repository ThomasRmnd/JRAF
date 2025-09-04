#include "selection/Muon.hpp"

MuonSelection::MuonSelection(const track& trk) :
    m_trk{trk}
{}

BasicMuonVetoSelection::BasicMuonVetoSelection(const track& trk, double radius, const TimeStamp& its, const TimeStamp& fts) :
    MuonSelection{trk},
    m_cyl{trk.ipos, trk.fpos, radius},
    m_win{trk.ts, its, fts}
{}

bool BasicMuonVetoSelection::isIn(const vertex& vtx) const {
    return m_cyl.isIn(vtx) && m_win.isIn(vtx);
}

EvolutiveCylindricalMuonVetoSelection::EvolutiveCylindricalMuonVetoSelection(const track& trk, double iradius, double tcoef) :
    MuonSelection{trk},
    m_iradius{iradius},
    m_tcoef{tcoef}
{}

bool EvolutiveCylindricalMuonVetoSelection::isIn(const vertex& vtx) const {
    if (vtx.ts < m_trk.ts) return false;
    double radius = m_iradius - static_cast<double>(vtx.ts - m_trk.ts) * m_tcoef;
    if (radius <= 0.0) {
        return false;
    }
    return CylinderVolumeSelection{m_trk.ipos, m_trk.fpos, radius}.isIn(vtx);
}

WaterPoolMuonVetoSelection::WaterPoolMuonVetoSelection(const track& trk, const TimeStamp& fts) :
    MuonSelection{trk},
    m_rts{trk.ts},
    m_fts{fts}
{}

bool WaterPoolMuonVetoSelection::isIn(const vertex& vtx) const {
    TimeStamp diff = vtx.ts - m_rts;
    return TimeStamp{0, 0} <= diff && diff <= m_fts;
}