#include "analysis/Analysis.hpp"

#include "SniperKernel/SniperLog.h"

Analysis::Analysis(const std::string& name, const std::string& method) :
    m_name{name},
    m_method{method},
    m_method_sel{method}
{}

bool Analysis::initialize() {
    m_tree = new TTree(m_name.c_str(), m_name.c_str());
    if (!m_tree) {
        std::cout << "[ERROR] Could not create the TTree for Analysis " << m_name << '\n';
        return false;
    }

    m_tree->Branch("posx_p", &posx_p);
    m_tree->Branch("posy_p", &posy_p);
    m_tree->Branch("posz_p", &posz_p);
    m_tree->Branch("e_p", &e_p);
    m_tree->Branch("totq_p", &totq_p);
    m_tree->Branch("sec_p", &sec_p);
    m_tree->Branch("nsec_p", &nsec_p);

    m_tree->Branch("posx_d", &posx_d);
    m_tree->Branch("posy_d", &posy_d);
    m_tree->Branch("posz_d", &posz_d);
    m_tree->Branch("e_d", &e_d);
    m_tree->Branch("totq_d", &totq_d);
    m_tree->Branch("sec_d", &sec_d);
    m_tree->Branch("nsec_d", &nsec_d);
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