#ifndef ANALYSISGROUPC_SELECTION_ENERGY_HPP_
#define ANALYSISGROUPC_SELECTION_ENERGY_HPP_

#include "selection/Selection.hpp"

class EnergySelection : public Selection{

public:

    virtual ~EnergySelection() = default;

    virtual bool isIn(const vertex& vtx) const = 0;

};

class IntervalEnergySelection : public EnergySelection {

public:

    IntervalEnergySelection(double lower, double upper);

    ~IntervalEnergySelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    double m_lower;
    double m_upper;

};

class MultiIntervalEnergySelection : public EnergySelection {

public:

    MultiIntervalEnergySelection(const std::vector<std::shared_ptr<IntervalEnergySelection>>& selecs);

    ~MultiIntervalEnergySelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    std::vector<std::shared_ptr<IntervalEnergySelection>> m_selecs;

};

#endif // ANALYSISGROUPC_SELECTION_ENERGY_HPP_