#ifndef ANALYSISGROUPC_EVENT_VERTEX_HPP_
#define ANALYSISGROUPC_EVENT_VERTEX_HPP_

#include <iostream>

#include "Context/TimeStamp.h"
#include "Event/RecVertex.h"

#include "utils/vec3.hpp"

struct vertex {

    vec3 pos;
    double energy;
    double totq;
    TimeStamp ts;
    std::string type;

    vertex(const vec3& pos_, double energy_, double totq_, const TimeStamp& ts_, const std::string& type_);

};

template<class _Char, class _Traits>
std::basic_ostream<_Char, _Traits>& operator<<(std::basic_ostream<_Char, _Traits>& os, const vertex& vtx) {
    return os << "pos: " << vtx.pos << ", energy: " << vtx.energy << ", totq: " << vtx.totq << ", ts: " << vtx.ts;
}

#endif // ANALYSISGROUPC_EVENT_VERTEX_HPP_