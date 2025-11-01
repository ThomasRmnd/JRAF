#include "event/Track.hpp"

#include <TVector3.h>

track::track(const std::string& method_, const JM::RecTrack& trk_, const TimeStamp& ts_, const loc& det_) :
    method{method_},
    ipos{trk_.start().x(), trk_.start().y(), trk_.start().z()},
    fpos{trk_.end().x(), trk_.end().y(), trk_.end().z()},
    totpe{static_cast<double>(trk_.peSum())},
    ts{ts_},
    det{det_},
    quality{static_cast<double>(trk_.quality())}
{}

track::track(const std::string& method_, const vec3& ipos_, const vec3& fpos_, double totpe_, const TimeStamp& ts_, const loc& det_, double quality_) :
    method{method_},
    ipos{ipos_},
    fpos{fpos_},
    totpe{totpe_},
    ts{ts_},
    det{det_},
    quality{quality_}
{}