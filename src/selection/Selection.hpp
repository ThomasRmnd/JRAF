#ifndef ANALYSISGROUPC_SELECTION_SELECTION_HPP_
#define ANALYSISGROUPC_SELECTION_SELECTION_HPP_

#include <functional>
#include <memory>
#include <vector>

#include "event/Vertex.hpp"

class Selection {

public:

    virtual ~Selection() = default;

    virtual bool isIn(const vertex& vtx) const = 0;

};

class MethodSelection : public Selection {

public:

    MethodSelection(const std::string& method);

    ~MethodSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    std::string m_method;

};

// ----------------------------------------------
// Logical classes 
// ----------------------------------------------

class LogicalAndSelection : public Selection {

public:

    LogicalAndSelection(const std::shared_ptr<Selection>& lhs, const std::shared_ptr<Selection>& rhs);

    ~LogicalAndSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    std::shared_ptr<Selection> m_lhs, m_rhs;

};

std::shared_ptr<Selection> logical_and(const std::shared_ptr<Selection>& lhs, const std::shared_ptr<Selection>& rhs);

class LogicalOrSelection : public Selection {

public:

    LogicalOrSelection(const std::shared_ptr<Selection>& lhs, const std::shared_ptr<Selection>& rhs);

    ~LogicalOrSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    std::shared_ptr<Selection> m_lhs, m_rhs;

};

std::shared_ptr<Selection> logical_or(const std::shared_ptr<Selection>& lhs, const std::shared_ptr<Selection>& rhs);

class LogicalNotSelection : public Selection {

public:

    LogicalNotSelection(const std::shared_ptr<Selection>& sel);

    ~LogicalNotSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    std::shared_ptr<Selection> m_sel;

};

std::shared_ptr<Selection> logical_not(const std::shared_ptr<Selection>& sel);

class FunctionalSelection : public Selection {

public:

    FunctionalSelection(const std::function<bool(const vertex&)>& func);

    ~FunctionalSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    std::function<bool(const vertex&)> m_func;

};

#endif // ANALYSISGROUPC_SELECTION_SELECTION_HPP_