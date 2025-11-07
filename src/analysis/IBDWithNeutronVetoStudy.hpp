#ifndef ANALYSISGROUPC_ANALYSIS_IBDWITHNEUTRONVETOSTUDY_HPP_
#define ANALYSISGROUPC_ANALYSIS_IBDWITHNEUTRONVETOSTUDY_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class IBDWithNeutronVetoStudy : public Analysis {

public:

    IBDWithNeutronVetoStudy(const std::string& name, const std::string& method);

    ~IBDWithNeutronVetoStudy() override = default;

    void process(const EventContext::View& events) override;

};

#endif // ANALYSISGROUPC_ANALYSIS_IBDWITHNEUTRONVETOSTUDY_HPP_