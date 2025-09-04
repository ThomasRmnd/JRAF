#ifndef ANALYSISGROUPC_ANALYSIS_FIRSTCROSSCHECKANALYSIS_HPP_
#define ANALYSISGROUPC_ANALYSIS_FIRSTCROSSCHECKANALYSIS_HPP_

#include "analysis/Analysis.hpp"

class FirstCrossCheckAnalysis : public Analysis {

public:

    FirstCrossCheckAnalysis(const std::string& name);

    ~FirstCrossCheckAnalysis() override = default;

    bool initialize() override;
    void process(JM::NavBuffer* buf) override;

private:

    double totq_p;
    double totq_d;

};

#endif // ANALYSISGROUPC_ANALYSIS_FIRSTCROSSCHECKANALYSIS_HPP_