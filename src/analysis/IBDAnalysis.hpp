#ifndef ANALYSISGROUPC_ANALYSIS_IBDANALYSIS_HPP_
#define ANALYSISGROUPC_ANALYSIS_IBDANALYSIS_HPP_

#include <vector>

#include "analysis/Analysis.hpp"
#include "event/IBD.hpp"

struct ibd_info {

    ibd pair;
    std::vector<vertex> neus;

    ibd_info(const vertex& prompt, const vertex& delayed) :
        pair{prompt, delayed}
    {}

};

class IBDAnalysis : public Analysis {

public:

    IBDAnalysis(const std::string& name, const std::string& method);

    ~IBDAnalysis() override = default;

    bool initialize() override;
    void process(const EventContext::View& events) override;

private:

    std::vector<double> posx_n;
    std::vector<double> posy_n;
    std::vector<double> posz_n;
    std::vector<double> e_n;
    std::vector<double> totq_n;
    std::vector<time_t> sec_n;
    std::vector<int> nsec_n;

    std::vector<std::string> method_mu;
    std::vector<int> loc_mu;
    std::vector<double> posx_mu;
    std::vector<double> posy_mu;
    std::vector<double> posz_mu;
    std::vector<double> dirx_mu;
    std::vector<double> diry_mu;
    std::vector<double> dirz_mu;
    std::vector<double> totq_mu;
    std::vector<time_t> sec_mu;
    std::vector<int> nsec_mu;
    std::vector<double> quality_mu;


};

#endif // ANALYSISGROUPC_ANALYSIS_IBDANALYSIS_HPP_