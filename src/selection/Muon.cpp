#include "selection/Muon.hpp"

MuonVetoSelection::MuonVetoSelection(const track& trk) :
    c_trk{trk}
{}

TimeRangeMuonVetoSelection::TimeRangeMuonVetoSelection(const track& trk, const TimeStamp& its, const TimeStamp& fts) :
    MuonVetoSelection{trk},
    m_trs{trk.ts, its, fts}
{}

bool TimeRangeMuonVetoSelection::isIn(const vertex& vtx) const {
    return m_trs.isIn(vtx);
}

CylindricalMuonVetoSelection::CylindricalMuonVetoSelection(const track& trk, double radius, const TimeStamp& its, const TimeStamp& fts) :
    MuonVetoSelection{trk},
    m_cyl{trk.ipos, trk.fpos, radius},
    m_trs{trk.ts, its, fts}
{}

bool CylindricalMuonVetoSelection::isIn(const vertex& vtx) const {
    return m_cyl.isIn(vtx) && m_trs.isIn(vtx);
}