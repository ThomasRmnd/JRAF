#include "selection/Volume.hpp"

FiducialVolumeSelection::FiducialVolumeSelection(double radius) :
    m_radius{radius}
{}

bool FiducialVolumeSelection::isIn(const vertex& vtx) const {
    return mag(vtx.pos) <= m_radius;
}

HeightVolumeSelection::HeightVolumeSelection(double lower_height, double upper_height) :
    m_lower_height{lower_height},
    m_upper_height{upper_height}
{}

bool HeightVolumeSelection::isIn(const vertex& vtx) const {
    return (m_lower_height <= vtx.pos.z) && (vtx.pos.z <= m_upper_height);
}

XYRadiusVolumeSelection::XYRadiusVolumeSelection(double radius_lower, double radius_upper) :
    m_radius_lower{radius_lower},
    m_radius_upper{radius_upper}
{}

bool XYRadiusVolumeSelection::isIn(const vertex& vtx) const {
    double r = std::sqrt(vtx.pos.x * vtx.pos.x + vtx.pos.y * vtx.pos.y);
    return (m_radius_lower <= r) && (r <= m_radius_upper);
}

SphereVolumeSelection::SphereVolumeSelection(const vec3& pos, double radius) :
    m_pos{pos},
    m_radius{radius}
{}

bool SphereVolumeSelection::isIn(const vertex& vtx) const {
    return mag(vtx.pos - m_pos) <= m_radius;
}

CylinderVolumeSelection::CylinderVolumeSelection(const vec3& ipos, const vec3& fpos, double radius) :
    m_ipos{ipos},
    m_fpos{fpos},
    m_radius{radius}
{}

bool CylinderVolumeSelection::isIn(const vertex& vtx) const {
    return mag(cross(unit(m_fpos - m_ipos), vtx.pos - m_ipos)) <= m_radius;
}