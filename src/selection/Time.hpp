#ifndef ANALYSISGROUPC_SELECTION_TIME_HPP_
#define ANALYSISGROUPC_SELECTION_TIME_HPP_

#include "selection/Selection.hpp"

#include "Context/TimeStamp.h"

class TimeSelection : public Selection {

public:

    virtual ~TimeSelection() = default;

    virtual bool isIn(const vertex& vtx) const = 0;

};

class UpperTimeSelection : public TimeSelection {

public:

    UpperTimeSelection(const TimeStamp& rts, const TimeStamp& uts);

    ~UpperTimeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    TimeStamp m_rts;
    TimeStamp m_uts;

};

class WindowTimeSelection : public TimeSelection {

public:

    WindowTimeSelection(const TimeStamp& rts, const TimeStamp& its, const TimeStamp& fts);

    ~WindowTimeSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    TimeStamp m_rts;
    TimeStamp m_its;
    TimeStamp m_fts;

};

#endif // ANALYSISGROUPC_SELECTION_TIME_HPP_