#ifndef ANALYSIS_EVENT_EVENTCONTEXT_HPP_
#define ANALYSIS_EVENT_EVENTCONTEXT_HPP_

#include <unordered_map>
#include <vector>

#include "EvtNavigator/NavBuffer.h"

#include "event/Track.hpp"
#include "event/Vertex.hpp"

class EventContextVertexView {

public:

    typedef std::vector<vertex>::const_iterator         const_iterator;
    typedef std::vector<vertex>::const_reverse_iterator const_reverse_iterator;

    EventContextVertexView(const_iterator first, const_iterator last) :
        m_first{first},
        m_last{last}
    {}

    const_iterator begin() const { return m_first; }
    const_iterator cbegin() const { return m_first; }
    const_reverse_iterator rbegin() const { return std::make_reverse_iterator(m_last); }
    const_reverse_iterator crbegin() const { return std::make_reverse_iterator(m_last); }

    const_iterator end() const { return m_last; }
    const_iterator cend() const { return m_last; }
    const_reverse_iterator rend() const { return std::make_reverse_iterator(m_first); }
    const_reverse_iterator crend() const { return std::make_reverse_iterator(m_first); }

private:

    const_iterator m_first, m_last;

};

class EventContextVertexRange {

public:

    typedef std::vector<vertex>::const_iterator         const_iterator;
    typedef std::vector<vertex>::const_reverse_iterator const_reverse_iterator;

    EventContextVertexRange(JM::NavBuffer* buf, const std::string& method);

    const_iterator begin() const { return m_vertices.begin(); }
    const_iterator cbegin() const { return m_vertices.cbegin(); }
    const_reverse_iterator rbegin() const { return m_vertices.rbegin(); }
    const_reverse_iterator crbegin() const { return m_vertices.crbegin(); }

    const_iterator end() const { return m_vertices.end(); }
    const_iterator cend() const { return m_vertices.cend(); }
    const_reverse_iterator rend() const { return m_vertices.rend(); }
    const_reverse_iterator crend() const { return m_vertices.crend(); }

    EventContextVertexView before() const { return EventContextVertexView{m_vertices.begin(), m_fcur}; }
    EventContextVertexView current() const { return EventContextVertexView{m_fcur, m_lcur}; }
    EventContextVertexView after() const { return EventContextVertexView{m_lcur, m_vertices.end()}; }
    const std::vector<vertex>& vertices() const { return m_vertices; }

private:

    std::vector<vertex> m_vertices;
    const_iterator m_fcur, m_lcur;

};

class EventContext {

public:

    class View {

    public:

        View(const std::vector<track>& tracks, const EventContextVertexRange& vertices) :
            c_tracks{tracks},
            c_vertices{vertices}
        {}

        EventContextVertexView before() const { return c_vertices.before(); }
        EventContextVertexView current() const { return c_vertices.current(); }
        EventContextVertexView after() const { return c_vertices.after(); }
        const std::vector<vertex>& vertices() const { return c_vertices.vertices(); }
        const std::vector<track>& tracks() const { return c_tracks; }

    private:

        const std::vector<track>& c_tracks;
        const EventContextVertexRange& c_vertices;

    };

    EventContext(JM::NavBuffer* buf, const std::vector<std::string>& methods);

    View view(const std::string& method) const { return View{m_tracks, m_vertices.at(method)}; }

private:

    std::vector<track> m_tracks;
    std::unordered_map<std::string, EventContextVertexRange> m_vertices;

};

#endif // ANALYSIS_EVENT_EVENTCONTEXT_HPP_