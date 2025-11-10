#ifndef MIXEDANALYSIS_ANALYSIS_ANALYSIS_HPP_
#define MIXEDANALYSIS_ANALYSIS_ANALYSIS_HPP_

#include "TFile.h"
#include "TTree.h"

#include "event/EventContext.hpp"
#include "selection/Selection.hpp"

class Analysis {

public:

    Analysis(const std::string& name, const std::string& method);

    virtual ~Analysis() = default;

    virtual bool initialize();
    const std::string& method() const;
    virtual void process(const EventContext::View& events) = 0;
    virtual bool write();

protected:

    std::string m_name;
    TTree* m_tree = nullptr;

    std::string m_method;
    MethodSelection m_method_sel;

    int run_id;

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

};

#endif // MIXEDANALYSIS_ANALYSIS_ANALYSIS_HPP_