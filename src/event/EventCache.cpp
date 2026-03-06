#include "event/EventCache.hpp"

#include "SniperKernel/SniperLog.h"

EventCache::CacheType EventCache::s_cache;
std::size_t EventCache::s_insert_counter = 0;
constexpr std::size_t EventCache::s_clean_interval;

std::shared_ptr<Event> EventCache::load(JM::EvtNavigator* nav) {
    TimeStamp ts{nav->TimeStamp().GetTimeSpec()};
    CacheType::iterator it = s_cache.find(ts);
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
    s_cache[ts] = evt;
    ++s_insert_counter;

    if (s_insert_counter % s_clean_interval == 0 && !s_cache.empty()) {
        clean(ts, TimeStamp{20, 0});
    }
}

bool EventCache::contains(JM::EvtNavigator* nav) {
    CacheType::iterator it = s_cache.find(TimeStamp{nav->TimeStamp().GetTimeSpec()});
    return (it != s_cache.end());
}

std::shared_ptr<Event> EventCache::get(JM::EvtNavigator* nav) {
    CacheType::iterator it = s_cache.find(TimeStamp{nav->TimeStamp().GetTimeSpec()});
    if (it != s_cache.end()) {
        return it->second;
    }
    return nullptr;
}

std::size_t EventCache::size() {
    return s_cache.size();
}

void EventCache::clean(const TimeStamp& newest_ts, const TimeStamp& window) {
    const TimeStamp allowed = newest_ts - window;

    EventCache::CacheType::iterator it = s_cache.lower_bound(allowed);
    if (it == s_cache.begin()) return;

    std::size_t removed = std::distance(s_cache.begin(), it);
    s_cache.erase(s_cache.begin(), it);

    std::cout << "EventCache: cleaned " << removed << " old events (older than "
              << window.GetSec() << "s, new size=" << s_cache.size() << ")\n";
}