#include "event/Track.hpp"

#include <TVector3.h>

track::track(const std::string& method_, const JM::RecTrack& trk_, const TimeStamp& ts_, const loc& det_) :
    method{method_},
    ipos{trk_.start().x(), trk_.start().y(), trk_.start().z()},
    fpos{trk_.end().x(), trk_.end().y(), trk_.end().z()},
    totq_cd{0.0},
    totq_wp{0.0},
    ts{ts_},
    det{det_},
    quality{static_cast<double>(trk_.quality())}
{
    if (det == loc::wp) {
        totq_wp = static_cast<double>(trk_.peSum());
    }
    else { // if det == loc::cd or other cases
        totq_cd = static_cast<double>(trk_.peSum());
    }
}

track::track(const std::string& method_, const vec3& ipos_, const vec3& fpos_, double totq_cd_, double totq_wp_, const TimeStamp& ts_, const loc& det_, double quality_) :
    method{method_},
    ipos{ipos_},
    fpos{fpos_},
    totq_cd{totq_cd_},
    totq_wp{totq_wp_},
    ts{ts_},
    det{det_},
    quality{quality_}
{}