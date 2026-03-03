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

};

#endif // MIXEDANALYSIS_ANALYSIS_ANALYSIS_HPP_