#ifndef ANALYSISGROUPC_SELECTION_TIME_HPP_
#define ANALYSISGROUPC_SELECTION_TIME_HPP_

#include "selection/Selection.hpp"

#include "Context/TimeStamp.h"

class TimeRangeSelection : public Selection {

public:

    TimeRangeSelection(const TimeStamp& rts, const TimeStamp& its, const TimeStamp& fts);

    ~TimeRangeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    TimeStamp m_rts, m_its, m_fts;

};

#endif // ANALYSISGROUPC_SELECTION_TIME_HPP_