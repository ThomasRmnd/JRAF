#include "event/Vertex.hpp"

vertex::vertex(
    const std::string& method_, const vec3& pos_, double energy_, const TimeStamp& ts_, 
    const calibration_context& calib_,
    const std::string& type_
) :
    method{method_},
    pos{pos_},
    energy{energy_},
    ts{ts_},
    calib{calib_},
    type{type_}
{}