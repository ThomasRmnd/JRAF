#ifndef ANALYSIS_EVENT_EVENTCACHE_HPP_
#define ANALYSIS_EVENT_EVENTCACHE_HPP_

#include <map>
#include <memory>
#include <unordered_map>

#include "EvtNavigator/EvtNavigator.h"

#include "event/Event.hpp"

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

private:
    // using CacheType = std::unordered_map<JM::EvtNavigator*, std::shared_ptr<Event>>;
    using CacheType = std::map<TimeStamp, std::shared_ptr<Event>>;

    static CacheType s_cache;

    // Private constructor to prevent instantiation.
    EventCache() = default;
};

#endif // ANALYSIS_EVENT_EVENTCACHE_HPP_
