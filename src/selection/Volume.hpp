#ifndef ANALYSISGROUPC_SELECTION_VOLUME_HPP_
#define ANALYSISGROUPC_SELECTION_VOLUME_HPP_

#include "UtilsThomas/math/vec3.hpp"

#include "selection/Selection.hpp"

class FiducialVolumeSelection : public Selection {

public:

    FiducialVolumeSelection(double thold);

    ~FiducialVolumeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    double m_thold2;

};

class HeigthRangeSelection : public Selection {

public:

    HeigthRangeSelection(double zmin, double zmax);

    ~HeigthRangeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    double m_zmin, m_zmax;

};

class RadialRangeSelection : public Selection {

public:

    RadialRangeSelection(double rhomin, double rhomax);

    ~RadialRangeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    double m_rhomin2, m_rhomax2;

};

class ChimneySelection : public Selection {

public:

    ChimneySelection(double z, double rho);

    ~ChimneySelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    HeigthRangeSelection m_bot;
    HeigthRangeSelection m_up;
    RadialRangeSelection m_rad;

};

class SphereSelection : public Selection {

public:

    SphereSelection(const vec3& pos, double radius);

    ~SphereSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    vec3 m_pos;
    double m_radius2;

};

class CylindricalSelection : public Selection {

public:

    CylindricalSelection(const vec3& ipos, const vec3& fpos, double radius);

    ~CylindricalSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    vec3 m_ipos, m_dir;
    double m_radius2;

};

#endif // ANALYSISGROUPC_SELECTION_VOLUME_HPP_