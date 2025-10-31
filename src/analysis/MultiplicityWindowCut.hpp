#ifndef ANALYSISGROUPC_ANALYSIS_MULTIPLICITYWINDOWCUT_HPP_
#define ANALYSISGROUPC_ANALYSIS_MULTIPLICITYWINDOWCUT_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class MultiplicityWindowCut : public Analysis {

public:

    MultiplicityWindowCut(const std::string& name, const std::string& method);

    ~MultiplicityWindowCut() override = default;

    bool initialize() override;
    void process(JM::NavBuffer* buf) override;

private:

    unsigned char m_window_type; // 0: before, 1: between, 2: after
    double posx_m;
    double posy_m;
    double posz_m;
    double e_m;
    double totq_m;
    time_t sec_m;
    int nsec_m;

};

#endif // ANALYSISGROUPC_ANALYSIS_MULTIPLICITYWINDOWCUT_HPP_