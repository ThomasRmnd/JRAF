#ifndef ANALYSISGROUPC_SELECTION_VERTEX_HPP_
#define ANALYSISGROUPC_SELECTION_VERTEX_HPP_

#include "selection/Selection.hpp"

#include "selection/Energy.hpp"
#include "selection/Muon.hpp"
#include "selection/Volume.hpp"
#include "selection/Time.hpp"

class VertexSelection : public Selection {

public:

    VertexSelection(const vertex& vtx);

    virtual ~VertexSelection() = default;

    virtual bool isIn(const vertex& vtx) const = 0;

    const vertex c_vtx;

};

class VertexCorrelationSelection : public VertexSelection {

public:

    VertexCorrelationSelection(const vertex& vtx, double radius, const TimeStamp& its, const TimeStamp& fts);

    ~VertexCorrelationSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    SphereSelection m_sph;
    TimeRangeSelection m_trs;

};

#endif // ANALYSISGROUPC_SELECTION_VERTEX_HPP_