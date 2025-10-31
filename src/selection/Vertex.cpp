#include "selection/Vertex.hpp"

VertexSelection::VertexSelection(const vertex& vtx) :
    c_vtx{vtx}
{}

VertexCorrelationSelection::VertexCorrelationSelection(const vertex& vtx, double radius, const TimeStamp& its, const TimeStamp& fts) :
    VertexSelection{vtx},
    m_sph{vtx.pos, radius},
    m_trs{vtx.ts, its, fts}
{}

bool VertexCorrelationSelection::isIn(const vertex& vtx) const {
    return m_sph.isIn(vtx) && m_trs.isIn(vtx);
}