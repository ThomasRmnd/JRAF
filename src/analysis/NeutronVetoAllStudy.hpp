#ifndef ANALYSISGROUPC_ANALYSIS_NEUTRONVETOALLSTUDY_HPP_
#define ANALYSISGROUPC_ANALYSIS_NEUTRONVETOALLSTUDY_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class NeutronVetoAllStudy : public Analysis {

public:

    NeutronVetoAllStudy(const std::string& name);

    ~NeutronVetoAllStudy() override = default;

    bool initialize() override;
    void process(JM::NavBuffer* buf) override;

private:

    double totq_p;
    double totq_d;

    double posx_e, posy_e, posz_e;
    double e_e;
    time_t sec_e;
    int nsec_e;
    double totq_e;

};

#endif // ANALYSISGROUPC_ANALYSIS_NEUTRONVETOALLSTUDY_HPP_