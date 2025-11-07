#ifndef ANALYSISGROUPC_ANALYSIS_FIRSTCROSSCHECKANALYSIS_HPP_
#define ANALYSISGROUPC_ANALYSIS_FIRSTCROSSCHECKANALYSIS_HPP_

#include "analysis/Analysis.hpp"

class FirstCrossCheckAnalysis : public Analysis {

public:

    FirstCrossCheckAnalysis(const std::string& name, const std::string& method);

    ~FirstCrossCheckAnalysis() override = default;

    void process(const EventContext::View& events) override;

};

#endif // ANALYSISGROUPC_ANALYSIS_FIRSTCROSSCHECKANALYSIS_HPP_