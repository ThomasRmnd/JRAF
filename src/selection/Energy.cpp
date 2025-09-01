#include "selection/Energy.hpp"

IntervalEnergySelection::IntervalEnergySelection(double lower, double upper) :
    m_lower{lower},
    m_upper{upper}
{}

bool IntervalEnergySelection::isIn(const vertex& vtx) const {
    return m_lower <= vtx.energy && vtx.energy <= m_upper;
}

MultiIntervalEnergySelection::MultiIntervalEnergySelection(const std::vector<std::shared_ptr<IntervalEnergySelection>>& selecs) :
    m_selecs(selecs)
{}

bool MultiIntervalEnergySelection::isIn(const vertex& vtx) const {
    for (const std::shared_ptr<IntervalEnergySelection>& selec : m_selecs) {
        if (selec->isIn(vtx)) return true;
    }
    return false;
}