#ifndef ANALYSISGROUPC_ANALYSIS_IBDWITHCYLINDRICALCUT_HPP_
#define ANALYSISGROUPC_ANALYSIS_IBDWITHCYLINDRICALCUT_HPP_

#include "analysis/Analysis.hpp"

#include "Context/TimeStamp.h"

class IBDWithCylindricalCut : public Analysis {

public:

    IBDWithCylindricalCut(const std::string& name, const std::string& method, double cyl_radius);

    ~IBDWithCylindricalCut() override = default;
    
    void process(JM::NavBuffer* buf) override;

private:

    double m_cyl_radius;

};

#endif // ANALYSISGROUPC_ANALYSIS_IBDWITHCYLINDRICALCUT_HPP_