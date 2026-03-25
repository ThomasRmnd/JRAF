#include "event/EventCache.hpp"

#include "SniperKernel/SniperLog.h"

EventCache::CacheType EventCache::s_cache;
#ifdef __EVENTCACHE_HASH_MAP_LONG64__
EventCache::OrderType EventCache::s_order;
#endif // __EVENTCACHE_HASH_MAP_LONG64__
std::size_t EventCache::s_insert_counter = 0;
constexpr std::size_t EventCache::s_clean_interval;

std::shared_ptr<Event> EventCache::load(JM::EvtNavigator* nav) {
    TimeStamp ts{nav->TimeStamp().GetTimeSpec()};
#ifndef __EVENTCACHE_HASH_MAP_LONG64__
    CacheType::iterator it = s_cache.find(ts);
#else
    std::int64_t key = ts_key(ts);
    CacheType::iterator it = s_cache.find(key);
#endif // __EVENTCACHE_HASH_MAP_LONG64__

    if (it != s_cache.end()) {
        return it->second;
    }

    std::shared_ptr<Event> evt = std::make_shared<Event>();
    if (!evt->load(nav)) {
        LogError << "EventCache::load: failed to load event for nav=" << nav << '\n';
        return nullptr;
    }
    insert(ts, evt);

    std::cout << "EventCache: loaded and cached event for nav=" << nav << '\n';
    return evt;
}

void EventCache::insert(const TimeStamp& ts, const std::shared_ptr<Event>& evt) {
#ifndef __EVENTCACHE_HASH_MAP_LONG64__
    s_cache[ts] = evt;
#else
    std::int64_t key = ts_key(ts);
    s_cache[key] = evt;
    s_order[key] = {};
#endif // __EVENTCACHE_HASH_MAP_LONG64__
    ++s_insert_counter;

    if (s_insert_counter % s_clean_interval == 0 && !s_cache.empty()) {
        clean(ts, TimeStamp{20, 0});
    }
}

bool EventCache::contains(JM::EvtNavigator* nav) {
#ifndef __EVENTCACHE_HASH_MAP_LONG64__
    return (s_cache.find(TimeStamp{nav->TimeStamp().GetTimeSpec()}) != s_cache.end());
#else
    return (s_cache.find(nav_key(nav)) != s_cache.end());
#endif // __EVENTCACHE_HASH_MAP_LONG64__
}

std::shared_ptr<Event> EventCache::get(JM::EvtNavigator* nav) {
#ifndef __EVENTCACHE_HASH_MAP_LONG64__
    CacheType::iterator it = s_cache.find(TimeStamp{nav->TimeStamp().GetTimeSpec()});
#else
    CacheType::iterator it = s_cache.find(nav_key(nav));
#endif // __EVENTCACHE_HASH_MAP_LONG64__
    if (it != s_cache.end()) {
        return it->second;
    }
    return nullptr;
}

std::size_t EventCache::size() {
    return s_cache.size();
}

void EventCache::clean(const TimeStamp& newest_ts, const TimeStamp& window) {
#ifndef __EVENTCACHE_HASH_MAP_LONG64__
    const TimeStamp allowed = newest_ts - window;

    CacheType::iterator boundary = s_cache.lower_bound(allowed);
    if (boundary == s_cache.begin()) return;

    std::size_t removed = std::distance(s_cache.begin(), boundary);
    s_cache.erase(s_cache.begin(), boundary);
#else
    const std::int64_t allowed = ts_key(newest_ts - window);

    OrderType::iterator boundary = s_order.lower_bound(allowed);
    if (boundary == s_order.begin()) return;

    std::size_t removed = 0ul;
    for (OrderType::iterator it = s_order.begin(); it != boundary; ++it) {
        s_cache.erase(it->first);
        ++removed;
    }
    s_order.erase(s_order.begin(), boundary);
#endif // __EVENTCACHE_HASH_MAP_LONG64__

    std::cout << "EventCache: cleaned " << removed << " old events (older than "
              << window.GetSec() << "s, new size=" << s_cache.size() << ")\n";
}