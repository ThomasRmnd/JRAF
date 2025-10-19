#ifndef ANALYSISGROUPC_ANALYSIS_CDWPALLCOSMOSTUDY_HPP_
#define ANALYSISGROUPC_ANALYSIS_CDWPALLCOSMOSTUDY_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class CdWpAllCosmoStudy : public Analysis {

public:

    CdWpAllCosmoStudy(const std::string& name);

    ~CdWpAllCosmoStudy() override = default;

    bool initialize() override;
    void process(JM::NavBuffer* buf) override;

private:

    double totq_p;
    double totq_d;

    double dlat_p, dlat_d;
    time_t dt_mu2p_sec, dt_mu2d_sec;
    int dt_mu2p_nsec, dt_mu2d_nsec;
};

#endif // ANALYSISGROUPC_ANALYSIS_CDWPALLCOSMOSTUDY_HPP_