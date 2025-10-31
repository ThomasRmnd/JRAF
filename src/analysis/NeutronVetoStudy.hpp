#ifndef ANALYSISGROUPC_ANALYSIS_NEUTRONVETOSTUDY_HPP_
#define ANALYSISGROUPC_ANALYSIS_NEUTRONVETOSTUDY_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class NeutronVetoStudy : public Analysis {

public:

    NeutronVetoStudy(const std::string& name, const std::string& method, double sph_radius, const TimeStamp& ts_window);

    ~NeutronVetoStudy() override = default;

    bool initialize() override;
    void process(JM::NavBuffer* buf) override;

private:

    double m_sph_radius;
    TimeStamp m_ts_window;

    double posx_e, posy_e, posz_e;
    double e_e;
    time_t sec_e;
    int nsec_e;
    double totq_e;

};

#endif // ANALYSISGROUPC_ANALYSIS_NEUTRONVETOSTUDY_HPP_