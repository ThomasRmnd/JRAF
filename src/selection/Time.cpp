#include "selection/Time.hpp"

TimeRangeSelection::TimeRangeSelection(const TimeStamp& rts, const TimeStamp& its, const TimeStamp& fts) :
    m_rts{rts},
    m_its{its},
    m_fts{fts}
{}

bool TimeRangeSelection::isIn(const vertex& vtx) const {
    TimeStamp diff = vtx.ts - m_rts;
    return m_its <= diff && diff <= m_fts;
}