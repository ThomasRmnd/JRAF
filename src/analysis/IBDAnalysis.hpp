#ifndef ANALYSISGROUPC_ANALYSIS_IBDANALYSIS_HPP_
#define ANALYSISGROUPC_ANALYSIS_IBDANALYSIS_HPP_

#include <vector>

#include "analysis/Analysis.hpp"
#include "event/IBD.hpp"

struct mult_info {
    vertex vtx;
    int type;
};

struct ibd_info {

    ibd pair;
    std::vector<vertex> neus;
    std::vector<mult_info> mults;

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
    std::vector<time_t> sec_n;
    std::vector<int> nsec_n;

    std::vector<double> totq_n;
    std::vector<double> meanq_n;
    std::vector<double> stdq_n;
    std::vector<double> minq_n;
    std::vector<double> maxq_n;
    std::vector<double> meant_n;
    std::vector<double> stdt_n;
    std::vector<std::size_t> npmt_n;
    std::vector<std::size_t> nhit_n;
    std::vector<double> meanhit_n;
    std::vector<double> stdhit_n;

    std::vector<double> posx_mult;
    std::vector<double> posy_mult;
    std::vector<double> posz_mult;
    std::vector<double> e_mult;
    std::vector<time_t> sec_mult;
    std::vector<int> nsec_mult;
    std::vector<int> mult_type; // 0 = before prompt, 1 = between, 2 = after delayed

    std::vector<double> totq_mult;
    std::vector<double> meanq_mult;
    std::vector<double> stdq_mult;
    std::vector<double> minq_mult;
    std::vector<double> maxq_mult;
    std::vector<double> meant_mult;
    std::vector<double> stdt_mult;
    std::vector<std::size_t> npmt_mult;
    std::vector<std::size_t> nhit_mult;
    std::vector<double> meanhit_mult;
    std::vector<double> stdhit_mult;

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