#include "selection/Selection.hpp"

CombinedSelection::CombinedSelection(const std::vector<std::shared_ptr<Selection>>& selections) :
    m_selections(selections)
{}

bool CombinedSelection::isIn(const vertex& vtx) const {
    for (const std::shared_ptr<Selection>& selection : m_selections) {
        if (!selection->isIn(vtx)) {
            return false;
        }
    }
    return true;
}