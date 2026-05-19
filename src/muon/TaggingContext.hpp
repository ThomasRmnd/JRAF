#ifndef JRAF_MUON_TAGGINGCONTEXT_HPP_
#define JRAF_MUON_TAGGINGCONTEXT_HPP_

#include "Context/TimeStamp.h"
#include "EvtNavigator/EvtNavigator.h"
#include "UtilsThomas/utils/DetectorType.hpp"

#include "event/CalibrationContext.hpp"

struct muon_tagging_context {

    TimeStamp ts;
    DetectorType det;

    calibration_context calid_cd;
    calibration_context calid_wp;

    JM::EvtNavigator* cd_nav = nullptr;
    JM::EvtNavigator* wp_nav = nullptr;

};

#endif // JRAF_MUON_TAGGINGCONTEXT_HPP_