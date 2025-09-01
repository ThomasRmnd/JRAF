#include "selection/Time.hpp"

UpperTimeSelection::UpperTimeSelection(const TimeStamp& rts, const TimeStamp& uts) :
    m_rts{rts},
    m_uts{uts}
{}

bool UpperTimeSelection::isIn(const vertex& vtx) const {
    TimeStamp diff = vtx.ts - m_rts;
    return TimeStamp{0, 0} <= diff && diff <= m_uts;
}

WindowTimeSelection::WindowTimeSelection(const TimeStamp& rts, const TimeStamp& its, const TimeStamp& fts) :
    m_rts{rts},
    m_its{its},
    m_fts{fts}
{}

bool WindowTimeSelection::isIn(const vertex& vtx) const {
    TimeStamp diff = vtx.ts - m_rts;
    return m_its <= diff && diff <= m_fts;
}