#ifndef ANALYSISGROUPC_SELECTION_MUON_HPP_
#define ANALYSISGROUPC_SELECTION_MUON_HPP_

#include "selection/Selection.hpp"

#include "event/Event.hpp"
#include "selection/Time.hpp"
#include "selection/Volume.hpp"

class MuonSelection : public Selection {

public:

    virtual ~MuonSelection() override = default;

    virtual bool isIn(const vertex& vtx) const = 0;

    const track* trkptr = nullptr;

};

class BasicMuonVetoSelection : public MuonSelection {

public:

    BasicMuonVetoSelection(const track& trk, double radius, const TimeStamp& its, const TimeStamp& fts);

    ~BasicMuonVetoSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    CylinderVolumeSelection m_cyl;
    WindowTimeSelection m_win;

};

class EvolutiveCylindricalMuonVetoSelection : public MuonSelection {

public:

    EvolutiveCylindricalMuonVetoSelection(const track& trk, double iradius, double tcoef);

    ~EvolutiveCylindricalMuonVetoSelection() override = default;
    
    bool isIn(const vertex& vtx) const override;

private:

    track m_trk;
    double m_iradius;
    double m_tcoef;

};

class WaterPoolMuonVetoSelection : public MuonSelection {

public:

    WaterPoolMuonVetoSelection(const track& trk, const TimeStamp& fts);

    ~WaterPoolMuonVetoSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    TimeStamp m_rts, m_fts;

};

#endif // ANALYSISGROUPC_SELECTION_MUON_HPP_