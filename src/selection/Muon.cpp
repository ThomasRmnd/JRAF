#include "selection/Muon.hpp"

BasicMuonVetoSelection::BasicMuonVetoSelection(const track& trk, double radius, const TimeStamp& its, const TimeStamp& fts) :
    m_cyl{trk.ipos, trk.fpos, radius},
    m_win{trk.ts, its, fts}
{
    trkptr = &trk;
}

bool BasicMuonVetoSelection::isIn(const vertex& vtx) const {
    return m_cyl.isIn(vtx) && m_win.isIn(vtx);
}

EvolutiveCylindricalMuonVetoSelection::EvolutiveCylindricalMuonVetoSelection(const track& trk, double iradius, double tcoef) :
    m_trk{trk},
    m_iradius{iradius},
    m_tcoef{tcoef}
{
    trkptr = &trk;
}

bool EvolutiveCylindricalMuonVetoSelection::isIn(const vertex& vtx) const {
    if (vtx.ts < m_trk.ts) return false;
    double radius = m_iradius - static_cast<double>(vtx.ts - m_trk.ts) * m_tcoef;
    if (radius <= 0.0) {
        return false;
    }
    return CylinderVolumeSelection{m_trk.ipos, m_trk.fpos, radius}.isIn(vtx);
}

WaterPoolMuonVetoSelection::WaterPoolMuonVetoSelection(const track& trk, const TimeStamp& fts) :
    m_rts{trk.ts},
    m_fts{fts}
{
    trkptr = &trk;
}

bool WaterPoolMuonVetoSelection::isIn(const vertex& vtx) const {
    TimeStamp diff = vtx.ts - m_rts;
    return TimeStamp{0, 0} <= diff && diff <= m_fts;
}