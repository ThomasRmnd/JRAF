#include "event/EventContext.hpp"

#include "event/EventCache.hpp"

EventContextVertexRange::EventContextVertexRange(JM::NavBuffer* buf, const std::string& method) {
    m_vertices.reserve(buf->size());
    std::ptrdiff_t fcur = 0l, lcur = 0l;
    for (JM::NavBuffer::Iterator it = buf->begin(); it != buf->end(); ++it) {
        JM::EvtNavigator* nav = it->get();
        std::shared_ptr<Event> evt = EventCache::load(nav);
        if (!evt) continue;
        if (it == buf->current()) {
            fcur = static_cast<std::ptrdiff_t>(m_vertices.size()); 
        }
        for (const vertex& vtx : evt->vertices) {
            if (vtx.method != method) continue;
            m_vertices.push_back(vtx);
        }
        if (it == buf->current()) {
            lcur = static_cast<std::ptrdiff_t>(m_vertices.size());
        }
    }
    m_fcur = m_vertices.begin() + fcur;
    m_lcur = m_vertices.begin() + lcur;
}

EventContext::EventContext(JM::NavBuffer* buf, const std::vector<std::string>& methods) {
    m_tracks.reserve(buf->size());
    for (JM::NavBuffer::Iterator it = buf->begin(); it != buf->end(); ++it) {
        JM::EvtNavigator* nav = it->get();
        std::shared_ptr<Event> evt = EventCache::load(nav);
        if (!evt) continue;
        m_tracks.insert(m_tracks.end(), evt->tracks.begin(), evt->tracks.end());
    }
    for (const std::string& method : methods) {
        m_vertices.emplace(method, EventContextVertexRange{buf, method});
    }
}