#include "selection/Energy.hpp"

EnergyRangeSelection::EnergyRangeSelection(double emin, double emax) :
    m_emin{emin},
    m_emax{emax}
{}

bool EnergyRangeSelection::isIn(const vertex& vtx) const {
    return m_emin <= vtx.energy && vtx.energy <= m_emax;
}

ChargeRangeSelection::ChargeRangeSelection(double qmin, double qmax) :
    m_qmin{qmin},
    m_qmax{qmax}
{}

bool ChargeRangeSelection::isIn(const vertex& vtx) const {
    return m_qmin <= vtx.totq && vtx.totq <= m_qmax;
}