#ifndef UTILS_EVENT_HPP_
#define UTILS_EVENT_HPP_

#include <iostream>

#include "utils/timestamp.hpp"
#include "utils/vec3.hpp"

struct vertex {

    vec3 pos;
    double e;
    timestamp ts;
    int run_id;

};

inline bool operator<(const vertex& lhs, const vertex& rhs) {
    return lhs.ts < rhs.ts;
}

inline std::ostream& operator<<(std::ostream& os, const vertex& v) {
    // return os << "Pos = " << v.pos << ", E = " << v.e << ", Q = " << v.q << ", Time = " << v.ts;
    return os << "E = " << v.e << ", Time = " << v.ts;
}

struct ibd {

    vertex prompt;
    vertex delayed;

};

inline bool operator<(const ibd& lhs, const ibd& rhs) {
    return lhs.prompt < rhs.prompt;
}

struct cosmogenic {

    vertex prompt;
    vertex delayed;
    double dlat_mu2p;
    double dlat_mu2d;
    double dt_mu2p;
    double dt_mu2d;

};

inline bool operator<(const cosmogenic& lhs, const cosmogenic& rhs) {
    return lhs.prompt < rhs.prompt;
}

struct vertex_metadata {

    double totq;
    double meanq;
    double stdq;
    double minq;
    double maxq;
    double meant;
    double stdt;
    std::size_t npmt;
    std::size_t nhit;
    double meanhit;
    double stdhit;

};

struct track {

    vec3 pos;
    vec3 dir;
    timestamp ts;
    double q;

};

#endif // UTILS_EVENT_HPP_