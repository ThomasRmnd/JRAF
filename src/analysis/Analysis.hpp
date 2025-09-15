#ifndef MIXEDANALYSIS_ANALYSIS_ANALYSIS_HPP_
#define MIXEDANALYSIS_ANALYSIS_ANALYSIS_HPP_

#include "TFile.h"
#include "TTree.h"

#include "EvtNavigator/NavBuffer.h"

class Analysis {

public:

    Analysis(const std::string& name);

    virtual ~Analysis() = default;

    virtual bool initialize();
    virtual void process(JM::NavBuffer* buf) = 0;
    virtual bool write();

protected:

    std::string m_name;
    TTree* m_tree = nullptr;

    double posx_p;
    double posy_p;
    double posz_p;
    double e_p;
    time_t sec_p;
    int nsec_p;

    double posx_d;
    double posy_d;
    double posz_d;
    double e_d;
    time_t sec_d;
    int nsec_d;

};

#endif // MIXEDANALYSIS_ANALYSIS_ANALYSIS_HPP_