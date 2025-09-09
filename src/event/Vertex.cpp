#include "event/Vertex.hpp"

vertex::vertex(const vec3& pos_, double energy_, double totq_, const TimeStamp& ts_, const std::string& type_) :
    pos{pos_},
    energy{energy_},
    totq{totq_},
    ts{ts_},
    type{type_}
{}