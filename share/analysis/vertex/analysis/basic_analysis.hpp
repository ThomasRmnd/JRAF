#ifndef ANALYSIS_BASIC_ANALYSIS_HPP_
#define ANALYSIS_BASIC_ANALYSIS_HPP_

#include "analysis/analysis.hpp"

class basic_analysis : public analysis_base {

public:

    basic_analysis(const std::string& name, const std::string& filepath, const std::string& suffix) :
        analysis_base{name}
    {
        std::string treename = "IBDAnalysis" + suffix;
        m_nav = navigator_manager::retrieve<basic_navigator>(filepath, treename);
        if (!m_nav->is_valid()) {
            std::cerr << "Cannot retrieve navigator of filepath " << filepath << " and treename " << treename << '\n';
            return;
        }
    }

    virtual ~basic_analysis() override = default;

    virtual std::shared_ptr<navigator_base> navigator() const override { return m_nav; }

    virtual bool selection() override = 0;
    virtual bool process() override = 0;

    virtual bool save() override = 0;
    virtual void result() override = 0;

protected:

    std::shared_ptr<basic_navigator> m_nav;

};

#endif // ANALYSIS_BASIC_ANALYSIS_HPP_