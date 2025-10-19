#ifndef ANALYSISGROUPC_ANALYSIS_IBDWITHNEUTRONVETOSTUDY_HPP_
#define ANALYSISGROUPC_ANALYSIS_IBDWITHNEUTRONVETOSTUDY_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class IBDWithNeutronVetoStudy : public Analysis {

public:

    IBDWithNeutronVetoStudy(const std::string& name);

    ~IBDWithNeutronVetoStudy() override = default;

    bool initialize() override;
    void process(JM::NavBuffer* buf) override;

private:

    double totq_p;
    double totq_d;

};

#endif // ANALYSISGROUPC_ANALYSIS_IBDWITHNEUTRONVETOSTUDY_HPP_