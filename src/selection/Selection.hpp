#ifndef ANALYSISGROUPC_SELECTION_SELECTION_HPP_
#define ANALYSISGROUPC_SELECTION_SELECTION_HPP_

#include <memory>
#include <vector>

#include "event/Vertex.hpp"

class Selection {

public:

    virtual ~Selection() = default;

    virtual bool isIn(const vertex& vtx) const = 0;

};

class CombinedSelection : public Selection {

public:

    CombinedSelection(const std::vector<std::shared_ptr<Selection>>& selections);

    ~CombinedSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    std::vector<std::shared_ptr<Selection>> m_selections;

};

#endif // ANALYSISGROUPC_SELECTION_SELECTION_HPP_