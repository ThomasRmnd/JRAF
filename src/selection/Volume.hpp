#ifndef ANALYSISGROUPC_SELECTION_VOLUME_HPP_
#define ANALYSISGROUPC_SELECTION_VOLUME_HPP_

#include "UtilsThomas/math/vec3.hpp"

#include "selection/Selection.hpp"

class VolumeSelection : public Selection {

public:

    virtual ~VolumeSelection() = default;

    virtual bool isIn(const vertex& vtx) const = 0;

};

class FiducialVolumeSelection : public VolumeSelection {

public:

    FiducialVolumeSelection(double radius);

    ~FiducialVolumeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    double m_radius;

};

class HeightVolumeSelection : public VolumeSelection {
public:

    HeightVolumeSelection(double lower_height, double upper_height);

    ~HeightVolumeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    double m_lower_height;
    double m_upper_height;

};

class XYRadiusVolumeSelection : public VolumeSelection {

public:

    XYRadiusVolumeSelection(double radius_lower, double radius_upper);

    ~XYRadiusVolumeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    double m_radius_lower;
    double m_radius_upper;

};

class SphereVolumeSelection : public VolumeSelection {

public:

    SphereVolumeSelection(const vec3& pos, double radius);
    
    ~SphereVolumeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    vec3 m_pos;
    double m_radius;

};

class CylinderVolumeSelection : public VolumeSelection {

public:

    CylinderVolumeSelection(const vec3& ipos, const vec3& fpos, double radius);

    ~CylinderVolumeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    vec3 m_ipos;
    vec3 m_fpos;
    double m_radius;

};

#endif // ANALYSISGROUPC_SELECTION_VOLUME_HPP_