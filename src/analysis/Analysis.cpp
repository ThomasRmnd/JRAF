#include "analysis/Analysis.hpp"

#include <iostream>

Analysis::Analysis(const std::string& name) :
    m_name(name)
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
    m_tree->Branch("sec_p", &sec_p);
    m_tree->Branch("nsec_p", &nsec_p);

    m_tree->Branch("posx_d", &posx_d);
    m_tree->Branch("posy_d", &posy_d);
    m_tree->Branch("posz_d", &posz_d);
    m_tree->Branch("e_d", &e_d);
    m_tree->Branch("sec_d", &sec_d);
    m_tree->Branch("nsec_d", &nsec_d);
    return true;
}

bool Analysis::write() {
    if (m_tree) {
        m_tree->Write();
    }
    return true;
}
