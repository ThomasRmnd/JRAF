#ifndef JRAF_EVENT_EVENTCACHE_HPP_
#define JRAF_EVENT_EVENTCACHE_HPP_

#include <map>
#include <memory>
#include <unordered_map>
#include <variant>

#include "EvtNavigator/EvtNavigator.h"

#include "event/Event.hpp"

#define __EVENTCACHE_HASH_MAP_LONG64__

/**
 * @class EventCache
 * @brief Global cache mapping JM::EvtNavigator* to loaded Event objects.
 *
 * This is a shared static cache that ensures each event is loaded only once.
 * Analyses can safely reuse events without redundant IO or reconstruction.
 */
class EventCache {

public:
    /**
     * @brief Get a shared_ptr to an Event for the given navigator.
     *
     * If the event is already cached, returns the existing one.
     * Otherwise, loads it and stores it in the cache.
     */
    static std::shared_ptr<Event> load(JM::EvtNavigator* nav);

    /**
     *  @brief Manually insert or update an event in the cache
     */
    static void insert(const TimeStamp& ts, const std::shared_ptr<Event>& evt);

    /**
     * @brief Check whether an event for the given navigator is already cached.
     */
    static bool contains(JM::EvtNavigator* nav);

    /**
     * @brief Get a shared_ptr to a cached Event without loading it.
     * Returns nullptr if not cached.
     */
    static std::shared_ptr<Event> get(JM::EvtNavigator* nav);

    /**
     * @brief For debugging: number of cached events.
     */
    static std::size_t size();

    /**
     * @brief Remove entries older than (ts - window)
     */
    static void clean(const TimeStamp& ts, const TimeStamp& window);

private:
#ifndef __EVENTCACHE_HASH_MAP_LONG64__
    using CacheType = std::map<TimeStamp, std::shared_ptr<Event>>;
    static CacheType s_cache;
#else
    using CacheType = std::unordered_map<std::int64_t, std::shared_ptr<Event>>;
    static CacheType s_cache;
    using OrderType = std::map<std::int64_t, std::monostate>;
    static OrderType s_order;
#endif // __EVENTCACHE_HASH_MAP_LONG64__

    static constexpr std::size_t s_clean_interval = 20000;
    static std::size_t s_insert_counter;

    // Private constructor to prevent instantiation.
    EventCache() = default;

#ifdef __EVENTCACHE_HASH_MAP_LONG64__
    static inline std::int64_t ts_key(const TimeStamp& ts) noexcept {
        return static_cast<std::int64_t>(ts.GetSec())  * 1'000'000'000LL + static_cast<std::int64_t>(ts.GetNanoSec());
    }

    static inline int64_t nav_key(JM::EvtNavigator* nav) noexcept {
        return ts_key(TimeStamp{nav->TimeStamp().GetTimeSpec()});
    }
#endif // __EVENTCACHE_HASH_MAP_LONG64__
};

#endif // JRAF_EVENT_EVENTCACHE_HPP_
