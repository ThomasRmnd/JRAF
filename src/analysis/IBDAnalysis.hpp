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
    bool write() override;

private:

    TTree* m_tree_cutflow = nullptr;

    int run_id;

    double posx_p;
    double posy_p;
    double posz_p;
    double e_p;
    time_t sec_p;
    int nsec_p;

    double totq_p;
    double meanq_p;
    double stdq_p;
    double minq_p;
    double maxq_p;
    double meant_p;
    double stdt_p;
    std::size_t npmt_p;
    std::size_t nhit_p;
    double meanhit_p;
    double stdhit_p;

    double posx_d;
    double posy_d;
    double posz_d;
    double e_d;
    time_t sec_d;
    int nsec_d;
    
    double totq_d;
    double meanq_d;
    double stdq_d;
    double minq_d;
    double maxq_d;
    double meant_d;
    double stdt_d;
    std::size_t npmt_d;
    std::size_t nhit_d;
    double meanhit_d;
    double stdhit_d;

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
    std::vector<double> iposx_mu;
    std::vector<double> iposy_mu;
    std::vector<double> iposz_mu;
    std::vector<double> fposx_mu;
    std::vector<double> fposy_mu;
    std::vector<double> fposz_mu;
    std::vector<double> totq_cd_mu;
    std::vector<double> totq_wp_mu;
    std::vector<time_t> sec_mu;
    std::vector<int> nsec_mu;
    std::vector<double> quality_mu;

    // Cut flow
    std::size_t cf_prompt_total = 0ul;
    std::size_t cf_prompt_fv = 0ul;
    std::size_t cf_prompt_energy = 0ul;
    std::size_t cf_prompt_muon = 0ul;
    std::size_t cf_pair_total = 0ul;
    std::size_t cf_pair_delayed_fv = 0ul;
    std::size_t cf_pair_delayed_energy = 0ul;
    std::size_t cf_pair_corr = 0ul;
    std::size_t cf_pair_delayed_muon = 0ul;
    std::size_t cf_ibd_final = 0ul;

};

#endif // ANALYSISGROUPC_ANALYSIS_IBDANALYSIS_HPP_