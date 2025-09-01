#include "event/Vertex.hpp"

vertex::vertex(const JM::RecVertex& vtx_, double totq_, const TimeStamp& ts_, const std::string& type_) :
    pos{vtx_.x(), vtx_.y(), vtx_.z()},
    energy{vtx_.energy()},
    totq{totq_},
    ts{ts_},
    type{type_}
{}