#ifndef ANALYSISGROUPC_EVENT_VERTEX_HPP_
#define ANALYSISGROUPC_EVENT_VERTEX_HPP_

#include <iostream>

#include "Context/TimeStamp.h"
#include "Event/RecVertex.h"
#include "UtilsThomas/math/vec3.hpp"

#include "event/CalibrationContext.hpp"

struct vertex {



    std::string method;
    vec3 pos;
    double energy;
    TimeStamp ts;

    calibration_context calib;

    std::string type;

    vertex(
        const std::string& method_, const vec3& pos_, double energy_, const TimeStamp& ts_, 
        const calibration_context& calib_, 
        const std::string& type_
    );

};

template<class _Char, class _Traits>
std::basic_ostream<_Char, _Traits>& operator<<(std::basic_ostream<_Char, _Traits>& os, const vertex& vtx) {
    return os << "method: " << vtx.method << ", pos: " << vtx.pos << ", energy: " << vtx.energy << ", ts: " << vtx.ts;
}

#endif // ANALYSISGROUPC_EVENT_VERTEX_HPP_