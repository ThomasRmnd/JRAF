#ifndef ANALYSISGROUPC_ANALYSIS_FIRSTCROSSCHECKANALYSIS_HPP_
#define ANALYSISGROUPC_ANALYSIS_FIRSTCROSSCHECKANALYSIS_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class FirstCrossCheckAnalysis : public Analysis {

public:

    FirstCrossCheckAnalysis(const std::string& name);

    ~FirstCrossCheckAnalysis() override = default;

    bool initialize() override;
    void process(JM::NavBuffer* buf) override;

private:

    double totq_p;
    double totq_d;

    TimeStamp m_prv_ts;
    time_t m_daq_sec;
    int m_daq_nsec;

    time_t m_muveto_sec;
    int m_muveto_nsec;

};

#endif // ANALYSISGROUPC_ANALYSIS_FIRSTCROSSCHECKANALYSIS_HPP_