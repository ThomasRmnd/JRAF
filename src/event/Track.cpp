#include "event/Track.hpp"

#include <TVector3.h>

track::track(const JM::RecTrack& trk_, const TimeStamp& ts_, const loc& det_) :
    ipos{trk_.start().x(), trk_.start().y(), trk_.start().z()},
    fpos{trk_.end().x(), trk_.end().y(), trk_.end().z()},
    totpe{trk_.peSum()},
    ts{ts_},
    det{det_}
{
    TVector3 ipos_(trk_.start().x(), trk_.start().y(), trk_.start().z());
    TVector3 fpos_(trk_.end().x(), trk_.end().y(), trk_.end().z());
    if ( std::abs( (fpos_ - ipos_).Unit().Cross(-ipos_).Mag() - mag(cross(unit(fpos - ipos), -ipos)) ) > 1.0 ) {
        std::cout << "[WARN] Clipping not matching between TVector3 (" << (fpos_ - ipos_).Unit().Cross(-ipos_).Mag() << ") and vec3 (" <<  mag(cross(unit(fpos - ipos), -ipos)) << '\n';
    }
}