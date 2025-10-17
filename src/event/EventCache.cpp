#include "event/EventCache.hpp"

#include "SniperKernel/SniperLog.h"

EventCache::CacheType EventCache::s_cache;

std::shared_ptr<Event> EventCache::load(JM::EvtNavigator* nav)
{
    if (!nav) {
        LogError << "EventCache::load: nullptr navigator" << std::endl;
        return nullptr;
    }

    CacheType::iterator it = s_cache.find(nav);
    if (it != s_cache.end()) {
        LogDebug << "EventCache: reusing cached event for nav=" << nav << std::endl;
        return it->second;
    }

    std::shared_ptr<Event> evt = std::make_shared<Event>();
    if (!evt->load(nav)) {
        LogError << "EventCache::load: failed to load event for nav=" << nav << std::endl;
        return nullptr;
    }

    s_cache[nav] = evt;
    LogDebug << "EventCache: loaded and cached event for nav=" << nav << std::endl;
    return evt;
}

bool EventCache::contains(JM::EvtNavigator* nav)
{
    CacheType::iterator it = s_cache.find(nav);
    return (it != s_cache.end());
}

std::shared_ptr<Event> EventCache::get(JM::EvtNavigator* nav)
{
    CacheType::iterator it = s_cache.find(nav);
    if (it != s_cache.end()) {
        return it->second;
    }
    return nullptr;
}

std::size_t EventCache::size()
{
    return s_cache.size();
}