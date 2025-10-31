#ifndef ANALYSISGROUPC_ANALYSIS_TTCOSMOSTUDY_HPP_
#define ANALYSISGROUPC_ANALYSIS_TTCOSMOSTUDY_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class TtCosmoStudy : public Analysis {

public:

    TtCosmoStudy(const std::string& name, const std::string& method, const TimeStamp& lwr_window, const TimeStamp& upr_window);

    ~TtCosmoStudy() override = default;

    bool initialize() override;
    void process(JM::NavBuffer* buf) override;

private:

    TimeStamp m_lwr_window, m_upr_window;

    double dlat_p, dlat_d;
    time_t dt_mu2p_sec, dt_mu2d_sec;
    int dt_mu2p_nsec, dt_mu2d_nsec;
};

#endif // ANALYSISGROUPC_ANALYSIS_TTCOSMOSTUDY_HPP_