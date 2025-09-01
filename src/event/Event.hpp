#ifndef ANALYSISGROUPC_EVENT_EVENT_HPP_
#define ANALYSISGROUPC_EVENT_EVENT_HPP_

#include <vector>

#include "Context/TimeStamp.h"
#include "EvtNavigator/EvtNavigator.h"

#include "event/Track.hpp"
#include "event/Vertex.hpp"

class Event {

public:

    Event() = default;

    bool load(JM::EvtNavigator* nav);

    std::vector<track> tracks;
    std::vector<vertex> vertices;
    double totq;
    TimeStamp ts;
    std::string type;

private:

    void loadCdTrack(JM::EvtNavigator* nav);
    void loadWpTrack(JM::EvtNavigator* nav);
    void loadCdVertex(JM::EvtNavigator* nav);
    void loadTrack(const JM::RecTrack* trk, const track::loc& det);
    void loadVertex(const JM::RecVertex* vtx);

};

template<class _Char, class _Traits>
std::basic_ostream<_Char, _Traits>& operator<<(std::basic_ostream<_Char, _Traits>& os, const Event& evt) {
    os << "Number of tracks: " << evt.tracks.size() << '\n';
    for (const track& trk : evt.tracks) {
        os << trk << '\n';
    }
    os << "Number of vertices: " << evt.vertices.size() << '\n';
    for (const vertex& vtx : evt.vertices) {
        os << vtx << '\n';
    }
    os << "timestamp: " << evt.ts;
    return os;
}

#endif // ANALYSISGROUPC_EVENT_EVENT_HPP_