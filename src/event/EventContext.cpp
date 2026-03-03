#include "event/EventContext.hpp"

#include "event/EventCache.hpp"

EventContext::EventContext(JM::NavBuffer* buf, const std::vector<std::string>& methods) {
    if (buf->curEvt()) m_run_id = buf->curEvt()->RunID();

    struct MethodData {
        std::vector<vertex> vertices;
        std::size_t fcur = 0;
        std::size_t lcur = 0;
    };

    std::unordered_map<std::string, MethodData> temp;
    temp.reserve(methods.size());
    for (const std::string& m : methods) {
        temp.emplace(m, MethodData{});
        temp[m].vertices.reserve(buf->size());
    }
    m_tracks.reserve(buf->size());

    for (JM::NavBuffer::Iterator it = buf->begin(); it != buf->end(); ++it) {
        JM::EvtNavigator* nav = it->get();
        std::shared_ptr<Event> evt = EventCache::load(nav);
        if (!evt) continue;
        m_tracks.insert(m_tracks.end(), evt->tracks.begin(), evt->tracks.end());

        bool is_current = (it == buf->current());

        for (const vertex& vtx : evt->vertices) {
            std::unordered_map<std::string, MethodData>::iterator found = temp.find(vtx.method);
            if (found == temp.end()) continue;

            MethodData& data = found->second;

            if (is_current) data.fcur = data.vertices.size();
            data.vertices.push_back(vtx);
            if (is_current) data.lcur = data.vertices.size();
        }
    }

    for (auto& [method, data] : temp) {
        m_vertices.emplace(
            method,
            EventContextVertexRange{
                std::move(data.vertices),
                data.fcur,
                data.lcur
            }
        );
    }
}