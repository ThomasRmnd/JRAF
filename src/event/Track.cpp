#include "event/Track.hpp"

#include <TVector3.h>

track::track(const std::string& method_, const JM::RecTrack& trk_, const TimeStamp& ts_, const loc& det_) :
    method{method_},
    ipos{trk_.start().x(), trk_.start().y(), trk_.start().z()},
    fpos{trk_.end().x(), trk_.end().y(), trk_.end().z()},
    totpe{trk_.peSum()},
    ts{ts_},
    det{det_},
    quality{trk_.quality()}
{}