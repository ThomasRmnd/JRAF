#ifndef ANALYSISGROUPC_ANALYSIS_NEUTRONANALYSIS_HPP_
#define ANALYSISGROUPC_ANALYSIS_NEUTRONANALYSIS_HPP_

#include "analysis/Analysis.hpp"

class NeutronAnalysis : public Analysis {

public:

    NeutronAnalysis(const std::string& name, const std::string& method);

    ~NeutronAnalysis() override = default;

    bool initialize() override;
    void process(const EventContext::View& events) override;

private:

    int run_id;

    double posx;
    double posy;
    double posz;
    double e;
    time_t sec;
    int nsec;

    double totq;
    double meanq;
    double stdq;
    double minq;
    double maxq;
    double meant;
    double stdt;
    std::size_t npmt;
    std::size_t nhit;
    double meanhit;
    double stdhit;

};

#endif // ANALYSISGROUPC_ANALYSIS_NEUTRONANALYSIS_HPP_