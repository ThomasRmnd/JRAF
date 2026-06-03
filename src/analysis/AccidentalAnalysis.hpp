#ifndef JRAF_ANALYSIS_ACCIDENTALANALYSIS_HPP_
#define JRAF_ANALYSIS_ACCIDENTALANALYSIS_HPP_

#include "analysis/Analysis.hpp"
#include "utils/event_info.hpp"

class AccidentalAnalysis : public Analysis {

public:

    AccidentalAnalysis(const std::string& name, const std::string& method);

    ~AccidentalAnalysis() override = default;

    bool initialize() override;
    void process(const EventContext::View& events) override;
    bool write() override;

private:

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

};

#endif // JRAF_ANALYSIS_ACCIDENTALANALYSIS_HPP_