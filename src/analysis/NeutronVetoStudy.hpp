#ifndef ANALYSISGROUPC_ANALYSIS_NEUTRONVETOSTUDY_HPP_
#define ANALYSISGROUPC_ANALYSIS_NEUTRONVETOSTUDY_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class NeutronVetoStudy : public Analysis {

public:

    NeutronVetoStudy(const std::string& name);

    ~NeutronVetoStudy() override = default;

    bool initialize() override;
    void process(JM::NavBuffer* buf) override;

private:

    double totq_p;
    double totq_d;

};

#endif // ANALYSISGROUPC_ANALYSIS_NEUTRONVETOSTUDY_HPP_