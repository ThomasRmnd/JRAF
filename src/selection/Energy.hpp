#ifndef JRAF_SELECTION_ENERGY_HPP_
#define JRAF_SELECTION_ENERGY_HPP_

#include "selection/Selection.hpp"

class EnergyRangeSelection : public Selection {

public:

    EnergyRangeSelection(double emin, double emax);

    ~EnergyRangeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    double m_emin, m_emax;

};

class ChargeRangeSelection : public Selection {

public:

    ChargeRangeSelection(double qmin, double qmax);

    ~ChargeRangeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    double m_qmin, m_qmax;

};

#endif // JRAF_SELECTION_ENERGY_HPP_