#ifndef UTILS_EVENT_HPP_
#define UTILS_EVENT_HPP_

#include <iostream>

#include "utils/timestamp.hpp"
#include "utils/vec3.hpp"

struct vertex {

    vec3 pos;
    double e;
    double q;
    timestamp ts;

};

std::ostream& operator<<(std::ostream& os, const vertex& v) {
    // os << "Pos = " << v.pos << ", E = " << v.e << ", Q = " << v.q << ", Time = " << v.ts;
    os << "E = " << v.e << ", Q = " << v.q << ", Time = " << v.ts;
    return os;
}

struct vertex_metadata {

    double meanq;
    double stdq;
    double minq;
    double maxq;
    std::size_t nhit;
    double meant;
    double stdt;

};

struct track {

    vec3 pos;
    vec3 dir;
    timestamp ts;
    double q;

};

#endif // UTILS_EVENT_HPP_