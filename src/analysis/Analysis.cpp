#include "analysis/Analysis.hpp"

#include "SniperKernel/SniperLog.h"

Analysis::Analysis(const std::string& name, const std::string& method) :
    m_name{name},
    m_method{method}
{}

bool Analysis::initialize() {
    m_tree = new TTree(m_name.c_str(), m_name.c_str());
    if (!m_tree) {
        std::cout << "[ERROR] Could not create the TTree for Analysis " << m_name << '\n';
        return false;
    }

    return true;
}

const std::string& Analysis::method() const {
    return m_method;
}

bool Analysis::write() {
    if (m_tree) {
        m_tree->Write();
    }
    return true;
}