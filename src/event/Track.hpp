#ifndef ANALYSISGROUPC_EVENT_TRACK_HPP_
#define ANALYSISGROUPC_EVENT_TRACK_HPP_

#include <iostream> 

#include "Context/TimeStamp.h"
#include "Event/RecTrack.h"
#include "UtilsThomas/math/vec3.hpp"

struct track {

    enum loc {
        none = 0,
        cd = 1,
        wp = 2,
        tt = 4
    };

    std::string method;
    vec3 ipos;
    vec3 fpos;
    double totpe;
    TimeStamp ts;
    loc det;
    float quality;

    track(const std::string& method_, const JM::RecTrack& trk_, const TimeStamp& ts_, const loc& det_);

};

template<class _Char, class _Traits>
std::basic_ostream<_Char, _Traits>& operator<<(std::basic_ostream<_Char, _Traits>& os, const track& trk) {
    return os << "method: " << method << ", ipos: " << trk.ipos << ", fpos: " << trk.fpos << ", totpe: " << trk.totpe << ", ts: " << trk.ts;
}

inline track::loc operator|(const track::loc& lhs, const track::loc& rhs) {
    return static_cast<track::loc>(static_cast<int>(lhs) | static_cast<int>(rhs));
}

inline track::loc operator&(const track::loc& lhs, const track::loc& rhs) {
    return static_cast<track::loc>(static_cast<int>(lhs) & static_cast<int>(rhs));
}

inline track::loc operator|=(track::loc& lhs, const track::loc& rhs) {
    lhs = lhs | rhs;
    return lhs;
}

inline track::loc operator&=(track::loc& lhs, const track::loc& rhs) {
    lhs = lhs & rhs;
    return lhs;
}

#endif // ANALYSISGROUPC_EVENT_TRACK_HPP_