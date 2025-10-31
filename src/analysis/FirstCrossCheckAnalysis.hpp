#ifndef ANALYSISGROUPC_ANALYSIS_FIRSTCROSSCHECKANALYSIS_HPP_
#define ANALYSISGROUPC_ANALYSIS_FIRSTCROSSCHECKANALYSIS_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class FirstCrossCheckAnalysis : public Analysis {

public:

    FirstCrossCheckAnalysis(const std::string& name, const std::string& method);

    ~FirstCrossCheckAnalysis() override = default;

    void process(JM::NavBuffer* buf) override;

};

#endif // ANALYSISGROUPC_ANALYSIS_FIRSTCROSSCHECKANALYSIS_HPP_