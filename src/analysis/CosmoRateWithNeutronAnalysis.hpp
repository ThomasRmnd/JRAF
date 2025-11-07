#ifndef ANALYSISGROUPC_ANALYSIS_COSMORATEWITHNEUTRONANALYSIS_HPP_
#define ANALYSISGROUPC_ANALYSIS_COSMORATEWITHNEUTRONANALYSIS_HPP_

#include "analysis/Analysis.hpp"

#include "event/Event.hpp"
#include "selection/Muon.hpp"

class CosmoRateWithNeutronAnalysis : public Analysis {

public:

    CosmoRateWithNeutronAnalysis(const std::string& name, const std::string& method);

    ~CosmoRateWithNeutronAnalysis() override = default;

    bool initialize() override;
    void process(const EventContext::View& events) override;

private:

    struct MuonAssociatedWithNeutron {

        TimeRangeMuonVetoSelection cosmo_veto;
        TimeRangeMuonVetoSelection neu_veto;
        std::vector<vertex> neu;

        MuonAssociatedWithNeutron(const TimeRangeMuonVetoSelection& cosmo_veto_, const TimeRangeMuonVetoSelection& neu_veto_) :
            cosmo_veto{cosmo_veto_},
            neu_veto{neu_veto_}
        {}

    };

    double dlat_p, dlat_d;
    time_t dt_mu2p_sec, dt_mu2d_sec;
    int dt_mu2p_nsec, dt_mu2d_nsec;

    double posx_n, posy_n, posz_n;
    double e_n;
    time_t sec_n;
    int nsec_n;
    double totq_n;

};

#endif // ANALYSISGROUPC_ANALYSIS_COSMORATEWITHNEUTRONANALYSIS_HPP_