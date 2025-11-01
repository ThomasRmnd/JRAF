#include "analysis/NavBufferCache.hpp"

#include "event/EventCache.hpp"

JM::NavBuffer* NavBufferCache::s_currentBuf = nullptr;
std::unordered_map<std::string, NavBufferCache::MethodData> NavBufferCache::s_data;

void NavBufferCache::prepare(JM::NavBuffer* buf, const std::vector<std::string>& methods) {
    if (s_currentBuf == buf && !s_data.empty()) return;

    clear();
    s_currentBuf = buf;

    for (const std::string& method : methods) {
        MethodData md;
        md.tracks.reserve(buf->size());
        md.vertices[VertexRegion::Before].reserve(buf->size() / 2);
        md.vertices[VertexRegion::After].reserve(buf->size() / 2);

        for (auto it = buf->begin(); it != buf->end(); ++it) {
            JM::EvtNavigator* nav = it->get();
            std::shared_ptr<Event> evt = EventCache::load(nav);
            if (!evt) continue;

            md.tracks.push_back(evt->tracks);

            for (const vertex& v : evt->vertices) {
                if (v.method != method) continue;
                if (it < buf->current())
                    md.vertices[VertexRegion::Before].push_back(v);
                else if (it > buf->current())
                    md.vertices[VertexRegion::After].push_back(v);
                else
                    md.vertices[VertexRegion::Current].push_back(v);
            }
        }

        s_data.emplace(method, std::move(md));
    }
}

const std::vector<vertex>& NavBufferCache::getVertices(const std::string& method, VertexRegion region) {
    static const std::vector<vertex> empty;
    std::unordered_map<std::string, NavBufferCache::MethodData>::iterator it = s_data.find(method);
    if (it == s_data.end()) return empty;
    std::map<NavBufferCache::VertexRegion, std::vector<vertex>>::iterator jt = it->second.vertices.find(region);
    if (jt == it->second.vertices.end()) return empty;
    return jt->second;
}

const std::vector<std::vector<track>>& NavBufferCache::getTracks(const std::string& method) {
    static const std::vector<std::vector<track>> empty;
    std::unordered_map<std::string, NavBufferCache::MethodData>::iterator it = s_data.find(method);
    if (it == s_data.end()) return empty;
    return it->second.tracks;
}

void NavBufferCache::clear() {
    s_data.clear();
    s_currentBuf = nullptr;
}
