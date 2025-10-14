#ifndef ANALYSISGROUPC_ANALYSIS_CDWPCOSMOSTUDY_HPP_
#define ANALYSISGROUPC_ANALYSIS_CDWPCOSMOSTUDY_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class CdWpCosmoStudy : public Analysis {

public:

    CdWpCosmoStudy(const std::string& name, const TimeStamp& lwr_window, const TimeStamp& upr_window);

    ~CdWpCosmoStudy() override = default;

    bool initialize() override;
    void process(JM::NavBuffer* buf) override;

private:

    TimeStamp m_lwr_window, m_upr_window;

    double totq_p;
    double totq_d;

    double dlat_p, dlat_d;
    time_t dt_mu2p_sec, dt_mu2d_sec;
    int dt_mu2p_nsec, dt_mu2d_nsec;
};

#endif // ANALYSISGROUPC_ANALYSIS_CDWPCOSMOSTUDY_HPP_