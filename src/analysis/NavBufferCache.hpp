#ifndef ANALYSIS_ANALYSIS_NAVBUFFERCACHE_HPP_
#define ANALYSIS_ANALYSIS_NAVBUFFERCACHE_HPP_

#include <map>
#include <vector>

#include "event/Event.hpp"

class NavBufferCache {
public:
    enum class VertexRegion { 
        Before, 
        Current, 
        After 
    };

    struct MethodData {
        std::vector<std::vector<track>> tracks;
        std::map<VertexRegion, std::vector<vertex>> vertices;
    };

    static void prepare(JM::NavBuffer* buf, const std::vector<std::string>& methods);

    static const std::vector<vertex>& getVertices(const std::string& method, VertexRegion region);

    static const std::vector<std::vector<track>>& getTracks(const std::string& method);

    static void clear();

private:

    static std::unordered_map<std::string, MethodData> s_data;
    static JM::NavBuffer* s_currentBuf;
};

#endif // ANALYSIS_ANALYSIS_NAVBUFFERCACHE_HPP_