#include "selection/Selection.hpp"

MethodSelection::MethodSelection(const std::string& method) :
    m_method{method}
{}

bool MethodSelection::isIn(const vertex& vtx) const {
    return vtx.method == m_method;
}

LogicalAndSelection::LogicalAndSelection(const std::shared_ptr<Selection>& lhs, const std::shared_ptr<Selection>& rhs) :
    m_lhs{lhs},
    m_rhs{rhs}
{}

bool LogicalAndSelection::isIn(const vertex& vtx) const {
    return m_lhs->isIn(vtx) && m_rhs->isIn(vtx);
}

std::shared_ptr<Selection> logical_and(const std::shared_ptr<Selection>& lhs, const std::shared_ptr<Selection>& rhs) {
    return std::make_shared<LogicalAndSelection>(lhs, rhs);
}

LogicalOrSelection::LogicalOrSelection(const std::shared_ptr<Selection>& lhs, const std::shared_ptr<Selection>& rhs) :
    m_lhs{lhs},
    m_rhs{rhs}
{}

bool LogicalOrSelection::isIn(const vertex& vtx) const {
    return m_lhs->isIn(vtx) || m_rhs->isIn(vtx);
}

std::shared_ptr<Selection> logical_or(const std::shared_ptr<Selection>& lhs, const std::shared_ptr<Selection>& rhs) {
    return std::make_shared<LogicalOrSelection>(lhs, rhs);
}

LogicalNotSelection::LogicalNotSelection(const std::shared_ptr<Selection>& sel) :
    m_sel{sel}
{}

bool LogicalNotSelection::isIn(const vertex& vtx) const {
    return !m_sel->isIn(vtx);
}

std::shared_ptr<Selection> logical_not(const std::shared_ptr<Selection>& sel) {
    return std::make_shared<LogicalNotSelection>(sel);
}

FunctionalSelection::FunctionalSelection(const std::function<bool(const vertex&)>& func) :
    m_func{func}
{}

bool FunctionalSelection::isIn(const vertex& vtx) const {
    return m_func(vtx);
}