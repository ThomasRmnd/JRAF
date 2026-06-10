#ifndef JRAF_VETO_VETO_HPP_
#define JRAF_VETO_VETO_HPP_

#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdTriggerHeader.h"
#include "Event/OecHeader.h"
#include "Event/WpCalibHeader.h"
#include "Event/WpTriggerHeader.h"
#include "EvtNavigator/EvtNavHelper.h"
#include "EvtNavigator/EvtNavigator.h"

#include "event/EventCache.hpp"

enum class VetoType : unsigned char {
    None,
    BeginningOfJob,
    MissingHeaders,
    BigGaps,
    Muon,
    MuonCd,
    MuonWp,
};

struct BeginningOfJobVetoTracker {

    bool check(std::size_t ievt) {
        return ievt == 1ul;
    }

};

struct MissingHeaderVetoTracker {

    JM::OecHeader* oec_hdr = nullptr;
    JM::CdLpmtCalibHeader* cd_lpmt_calib_hdr = nullptr;
    JM::CdTriggerHeader* cd_trig_hdr = nullptr;
    JM::WpCalibHeader* wp_calib_hdr = nullptr;
    JM::WpTriggerHeader* wp_trig_hdr = nullptr;

    bool check(JM::EvtNavigator* nav) {
        oec_hdr = JM::getHeaderObject<JM::OecHeader>(nav);
        cd_lpmt_calib_hdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
        cd_trig_hdr = JM::getHeaderObject<JM::CdTriggerHeader>(nav);
        wp_calib_hdr = JM::getHeaderObject<JM::WpCalibHeader>(nav);
        wp_trig_hdr = JM::getHeaderObject<JM::WpTriggerHeader>(nav);
        if (
            !oec_hdr ||
            (!cd_lpmt_calib_hdr && !wp_calib_hdr) ||
            (!cd_trig_hdr && !wp_trig_hdr)
        ) {
            return true;
        }
        return false;
    }

};

struct BigGapsVetoTracker {

    TimeStamp cd_thold{0, 50000000};
    TimeStamp wp_thold{0, 70000000};
    TimeStamp last_cd_ts{0, 0};
    TimeStamp last_wp_ts{0, 0};
    bool is_cd_initialized = false;
    bool is_wp_initialized = false;

    JM::CdLpmtCalibHeader* cd_lpmt_calib_hdr = nullptr;
    JM::WpCalibHeader* wp_calib_hdr = nullptr;

    bool check(JM::EvtNavigator* nav) {
        TimeStamp ts{nav->TimeStamp().GetTimeSpec()};
        
        cd_lpmt_calib_hdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
        wp_calib_hdr = JM::getHeaderObject<JM::WpCalibHeader>(nav);
        
        if (cd_lpmt_calib_hdr) {
            if (!is_cd_initialized) {
                last_cd_ts = ts;
                is_cd_initialized = true;
                return false;
            }
            TimeStamp diff = ts - last_cd_ts;
            last_cd_ts = ts;
            return diff > cd_thold;
        }
        else if (wp_calib_hdr) {
            if (!is_wp_initialized) {
                last_wp_ts = ts;
                is_wp_initialized = true;
                return false;
            }
            TimeStamp diff = ts - last_wp_ts;
            last_wp_ts = ts;
            return diff > wp_thold;
        }
        return false;
    }

};

struct MuonVetoTracker {

    bool check(JM::EvtNavigator* nav) {
        std::shared_ptr<Event> evt = EventCache::load(nav);
        return !evt->tracks.empty();
    }

    VetoType type(JM::EvtNavigator* nav) {
        std::shared_ptr<Event> evt = EventCache::load(nav);
        if (evt->tracks.empty()) {
            return VetoType::None;
        }
        if (evt->totq_cd > 0.0) {
            return VetoType::MuonCd;
        }
        if (evt->totq_wp > 0.0) {
            return VetoType::MuonWp;
        }
        return VetoType::Muon;
    }

};

#endif // JRAF_VETO_VETO_HPP_