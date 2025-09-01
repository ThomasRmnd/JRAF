#ifndef ANALYSISGROUPC_EVENT_TRACK_HPP_
#define ANALYSISGROUPC_EVENT_TRACK_HPP_

#include <iostream> 

#include "Context/TimeStamp.h"
#include "Event/RecTrack.h"

#include "utils/vec3.hpp"

struct track {

    enum loc {
        none = 0,
        cd = 1,
        wp = 2
    };

    vec3 ipos;
    vec3 fpos;
    double totpe;
    TimeStamp ts;
    loc det;

    track(const JM::RecTrack& trk_, const TimeStamp& ts_, const loc& det_);

};

template<class _Char, class _Traits>
std::basic_ostream<_Char, _Traits>& operator<<(std::basic_ostream<_Char, _Traits>& os, const track& trk) {
    return os << "ipos: " << trk.ipos << ", fpos: " << trk.fpos << ", totpe: " << trk.totpe << ", ts: " << trk.ts;
}

#endif // ANALYSISGROUPC_EVENT_TRACK_HPP_