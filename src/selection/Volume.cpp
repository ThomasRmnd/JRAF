#include "selection/Volume.hpp"

#include <limits>

FiducialVolumeSelection::FiducialVolumeSelection(double thold) :
    m_thold2{thold * thold}
{}

bool FiducialVolumeSelection::isIn(const vertex& vtx) const {
    return mag2(vtx.pos) <= m_thold2;
}

HeigthRangeSelection::HeigthRangeSelection(double zmin, double zmax) :
    m_zmin{zmin},
    m_zmax{zmax}
{}

bool HeigthRangeSelection::isIn(const vertex& vtx) const {
    return m_zmin <= vtx.pos.z && vtx.pos.z <= m_zmax;
}

RadialRangeSelection::RadialRangeSelection(double rhomin, double rhomax) :
    m_rhomin2{rhomin * rhomin},
    m_rhomax2{rhomax * rhomax}
{}

bool RadialRangeSelection::isIn(const vertex& vtx) const {
    double rho = vtx.pos.x * vtx.pos.x + vtx.pos.y * vtx.pos.y;
    return m_rhomin2 <= rho && rho <= m_rhomax2;
}

ChimneySelection::ChimneySelection(double z, double rho) :
    m_bot{-std::numeric_limits<double>::max(), -z},
    m_up{z, std::numeric_limits<double>::max()},
    m_rad{0.0, rho}
{}

bool ChimneySelection::isIn(const vertex& vtx) const {
    return (m_bot.isIn(vtx) && m_rad.isIn(vtx)) || (m_up.isIn(vtx) && m_rad.isIn(vtx));
}

SphereSelection::SphereSelection(const vec3& pos, double radius) :
    m_pos{pos},
    m_radius2{radius * radius}
{}

bool SphereSelection::isIn(const vertex& vtx) const {
    return mag2(vtx.pos - m_pos) <= m_radius2;
}

CylindricalSelection::CylindricalSelection(const vec3& ipos, const vec3& fpos, double radius) :
    m_ipos{ipos},
    m_dir{unit(fpos - ipos)},
    m_radius2{radius * radius}
{}

bool CylindricalSelection::isIn(const vertex& vtx) const {
    return mag2(cross(m_dir, vtx.pos - m_ipos)) <= m_radius2;
}