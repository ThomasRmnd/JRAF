#ifndef MIXEDANALYSIS_ANALYSIS_ANALYSIS_HPP_
#define MIXEDANALYSIS_ANALYSIS_ANALYSIS_HPP_

#include "TFile.h"
#include "TTree.h"

#include "EvtNavigator/NavBuffer.h"

#include "selection/Selection.hpp"

class Analysis {

public:

    Analysis(const std::string& name, const std::string& method);

    virtual ~Analysis() = default;

    virtual bool initialize();
    virtual void process(JM::NavBuffer* buf) = 0;
    virtual bool write();

protected:

    std::string m_name;
    TTree* m_tree = nullptr;

    MethodSelection m_method_sel;

    double posx_p;
    double posy_p;
    double posz_p;
    double e_p;
    double totq_p;
    time_t sec_p;
    int nsec_p;

    double posx_d;
    double posy_d;
    double posz_d;
    double e_d;
    double totq_d;
    time_t sec_d;
    int nsec_d;

    void extractEvent(JM::NavBuffer* buf, std::vector<std::vector<track>>& tracks, std::vector<vertex>& cur_vertices, std::vector<vertex>& bef_vertices, std::vector<vertex>& aft_vertices);

};

#endif // MIXEDANALYSIS_ANALYSIS_ANALYSIS_HPP_