#ifndef EVENT_HPP_
#define EVENT_HPP_

#include "timestamp.hpp"
#include "vec3.hpp"

struct vertex {

    vec3 pos;
    timestamp ts;
    double e;
    double q;

};

struct track {

    vec3 pos;
    vec3 dir;
    timestamp ts;
    double q;

};

#endif // EVENT_HPP_