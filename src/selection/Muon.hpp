#ifndef ANALYSISGROUPC_SELECTION_MUON_HPP_
#define ANALYSISGROUPC_SELECTION_MUON_HPP_

#include "selection/Selection.hpp"

#include "event/Event.hpp"
#include "selection/Time.hpp"
#include "selection/Volume.hpp"

class MuonVetoSelection : public Selection {

public:

    MuonVetoSelection(const track& trk);

    virtual ~MuonVetoSelection() = default;

    virtual bool isIn(const vertex& vtx) const = 0;

    const track c_trk;

};

class TimeRangeMuonVetoSelection : public MuonVetoSelection {

public:

    TimeRangeMuonVetoSelection(const track& trk, const TimeStamp& its, const TimeStamp& fts);

    ~TimeRangeMuonVetoSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    TimeRangeSelection m_trs;

};

class CylindricalMuonVetoSelection : public MuonVetoSelection {

public:

    CylindricalMuonVetoSelection(const track& trk, double radius, const TimeStamp& its, const TimeStamp& fts);

    ~CylindricalMuonVetoSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    CylindricalSelection m_cyl; 
    TimeRangeSelection m_trs;   

};

#endif // ANALYSISGROUPC_SELECTION_MUON_HPP_