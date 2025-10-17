#include "event/EventCache.hpp"

#include "SniperKernel/SniperLog.h"

EventCache::CacheType EventCache::s_cache;

std::shared_ptr<Event> EventCache::load(JM::EvtNavigator* nav)
{
    if (!nav) {
        LogError << "EventCache::load: nullptr navigator\n";
        return nullptr;
    }

    // Check if it exists in cache
    auto it = s_cache.find(nav);
    if (it != s_cache.end()) {
        if (auto existing = it->second.lock()) {
            LogInfo << "EventCache: reusing cached event for nav=" << nav << '\n';
            return existing;
        }
    }

    // Otherwise, load it
    auto evt = std::make_shared<Event>();
    if (!evt->load(nav)) {
        LogError << "EventCache::load: failed to load event for nav=" << nav << '\n';
        return nullptr;
    }

    s_cache[nav] = evt;
    LogInfo << "EventCache: loaded and cached event for nav=" << nav << '\n';
    return evt;
}

bool EventCache::contains(JM::EvtNavigator* nav)
{
    auto it = s_cache.find(nav);
    return (it != s_cache.end() && !it->second.expired());
}

std::shared_ptr<Event> EventCache::get(JM::EvtNavigator* nav)
{
    auto it = s_cache.find(nav);
    if (it != s_cache.end()) {
        return it->second.lock();
    }
    return nullptr;
}

std::size_t EventCache::size()
{
    return s_cache.size();
}