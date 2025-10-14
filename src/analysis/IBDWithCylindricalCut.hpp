#ifndef ANALYSISGROUPC_ANALYSIS_IBDWITHCYLINDRICALCUT_HPP_
#define ANALYSISGROUPC_ANALYSIS_IBDWITHCYLINDRICALCUT_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class IBDWithCylindricalCut : public Analysis {

public:

    IBDWithCylindricalCut(const std::string& name, double cyl_radius);

    ~IBDWithCylindricalCut() override = default;

    bool initialize() override;
    void process(JM::NavBuffer* buf) override;

private:

    double m_cyl_radius;

    double totq_p;
    double totq_d;

};

#endif // ANALYSISGROUPC_ANALYSIS_IBDWITHCYLINDRICALCUT_HPP_