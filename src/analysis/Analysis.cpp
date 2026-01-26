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
    m_tree_cutflow = new TTree((m_name + "__CutFlow").c_str(), (m_name + "__CutFlow").c_str());
    if (!m_tree_cutflow) {
        std::cout << "[ERROR] Could not create the TTree for CutFlow of Analysis " << m_name << '\n';
        return false;
    }

    m_tree->Branch("posx_p", &posx_p);
    m_tree->Branch("posy_p", &posy_p);
    m_tree->Branch("posz_p", &posz_p);
    m_tree->Branch("e_p", &e_p);
    m_tree->Branch("sec_p", &sec_p);
    m_tree->Branch("nsec_p", &nsec_p);

    m_tree->Branch("totq_p", &totq_p);
    m_tree->Branch("meanq_p", &meanq_p);
    m_tree->Branch("stdq_p", &stdq_p);
    m_tree->Branch("minq_p", &minq_p);
    m_tree->Branch("maxq_p", &maxq_p);
    m_tree->Branch("meant_p", &meant_p);
    m_tree->Branch("stdt_p", &stdt_p);
    m_tree->Branch("npmt_p", &npmt_p);
    m_tree->Branch("nhit_p", &nhit_p);
    m_tree->Branch("meanhit_p", &meanhit_p);
    m_tree->Branch("stdhit_p", &stdhit_p);

    m_tree->Branch("posx_d", &posx_d);
    m_tree->Branch("posy_d", &posy_d);
    m_tree->Branch("posz_d", &posz_d);
    m_tree->Branch("e_d", &e_d);
    m_tree->Branch("sec_d", &sec_d);
    m_tree->Branch("nsec_d", &nsec_d);

    m_tree->Branch("totq_d", &totq_d);
    m_tree->Branch("meanq_d", &meanq_d);
    m_tree->Branch("stdq_d", &stdq_d);
    m_tree->Branch("minq_d", &minq_d);
    m_tree->Branch("maxq_d", &maxq_d);
    m_tree->Branch("meant_d", &meant_d);
    m_tree->Branch("stdt_d", &stdt_d);
    m_tree->Branch("npmt_d", &npmt_d);
    m_tree->Branch("nhit_d", &nhit_d);
    m_tree->Branch("meanhit_d", &meanhit_d);
    m_tree->Branch("stdhit_d", &stdhit_d);

    return true;
}

const std::string& Analysis::method() const {
    return m_method;
}

bool Analysis::write() {
    if (m_tree) {
        m_tree->Write();
    }
    if (m_tree_cutflow) {
        m_tree_cutflow->Fill();
        m_tree_cutflow->Write();
    }
    return true;
}