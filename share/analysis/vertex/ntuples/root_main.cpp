#include <algorithm>
#include <chrono>
#include <iostream>
#include <set>
#include <string>

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include "analysis/basic_analysis.hpp"
#include "analysis/cosmo_rate_analysis.hpp"
#include "analysis/cosmo_rate_neutron_veto_analysis.hpp"
#include "analysis/cosmo_shape_muon_changing_veto_analysis.hpp"
#include "analysis/cosmo_shape_muon_standard_analysis.hpp"
#include "analysis/cosmo_shape_muon_with_neutron_analysis.hpp"
#include "analysis/cosmo_shape_neutron_analysis.hpp"
#include "analysis/ibd_no_neutron_veto_analysis.hpp"
#include "analysis/ibd_no_neutron_veto_muon_veto_analysis.hpp"
#include "analysis/ibd_standard_analysis.hpp"
#include "analysis/ibd_standard_muon_veto_analysis.hpp"
#include "analysis/ibd_standard_muon_with_neutron_veto_analysis.hpp"

int root_main(const std::string& filepath) {
    std::string suffix = "__OMILREC_JVtx";

    analysis_registry registry;
    analysis_manager manager(registry);



    std::shared_ptr<analysis_base> ibd_no_neutron_veto_analysis_omilrec_jvertex(
        new ibd_no_neutron_veto_analysis(
            "ibd_no_neutron_veto_analysis_omilrec_jvertex", 
            filepath, 
            suffix
        )
    );
    if (!registry.book(ibd_no_neutron_veto_analysis_omilrec_jvertex)) return 1;

    std::shared_ptr<analysis_base> ibd_no_neutron_veto_muon_veto_analysis_omilrec_jvertex(
        new ibd_no_neutron_veto_muon_veto_analysis(
            "ibd_no_neutron_veto_muon_veto_analysis_omilrec_jvertex", 
            filepath, 
            suffix, 
            "CdWpTtChi2", 
            timestamp{0, 5000000}, 
            timestamp{0, 1200000000}, 
            3000.0
        )
    );
    if (!registry.book(ibd_no_neutron_veto_muon_veto_analysis_omilrec_jvertex)) return 1;

    std::shared_ptr<analysis_base> ibd_standard_analysis_omilrec_jvertex(
        new ibd_standard_analysis(
            "ibd_standard_analysis_omilrec_jvertex", 
            filepath, 
            suffix
        )
    );
    if (!registry.book(ibd_standard_analysis_omilrec_jvertex)) return 1;

    std::shared_ptr<analysis_base> ibd_standard_muon_veto_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex(
        new ibd_standard_muon_veto_analysis(
            "ibd_standard_muon_veto_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex", 
            filepath, 
            suffix, 
            "CdWpTtChi2", 
            timestamp{0, 5000000}, 
            timestamp{0, 1200000000}, 
            3000.0
        )
    );
    if (!registry.book(ibd_standard_muon_veto_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex)) return 1;

    std::shared_ptr<analysis_base> ibd_standard_muon_with_neutron_veto_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex(
        new ibd_standard_muon_with_neutron_veto_analysis(
            "ibd_standard_muon_with_neutron_veto_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex", 
            filepath, 
            suffix, 
            "CdWpTtChi2", 
            timestamp{0, 5000000}, 
            timestamp{0, 1200000000}, 
            3000.0
        )
    );
    if (!registry.book(ibd_standard_muon_with_neutron_veto_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex)) return 1;

    std::shared_ptr<analysis_base> ibd_standard_muon_with_neutron_veto_analysis_cdwpttchi2_1m_0_5s_omilrec_jvertex(
        new ibd_standard_muon_with_neutron_veto_analysis(
            "ibd_standard_muon_with_neutron_veto_analysis_cdwpttchi2_1m_0_5s_omilrec_jvertex", 
            filepath, 
            suffix, 
            "CdWpTtChi2", 
            timestamp{0, 5000000}, 
            timestamp{0, 500000000}, 
            1000.0
        )
    );
    if (!registry.book(ibd_standard_muon_with_neutron_veto_analysis_cdwpttchi2_1m_0_5s_omilrec_jvertex)) return 1;



    // std::shared_ptr<analysis_base> cosmo_rate_analysis_omilrec_jvertex(new cosmo_rate_analysis("cosmo_rate_analysis_omilrec_jvertex", filepath, suffix));
    // if (!registry.book(cosmo_rate_analysis_omilrec_jvertex)) return 1;

    // std::shared_ptr<analysis_base> cosmo_rate_neutron_veto_analysis_omilrec_jvertex(new cosmo_rate_neutron_veto_analysis("cosmo_rate_neutron_veto_analysis_omilrec_jvertex", filepath, suffix));
    // if (!registry.book(cosmo_rate_neutron_veto_analysis_omilrec_jvertex)) return 1;



    std::shared_ptr<analysis_base> cosmo_shape_muon_standard_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex(
        new cosmo_shape_muon_standard_analysis(
            "cosmo_shape_muon_standard_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex", 
            filepath, 
            suffix, 
            "CdWpTtChi2", 
            timestamp{0, 5000000}, 
            timestamp{0, 1200000000}, 
            timestamp{0, -1200000000}, timestamp{0, -5000000}, 3000.0));
    if (!registry.book(cosmo_shape_muon_standard_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex)) return 1;

    std::shared_ptr<analysis_base> cosmo_shape_muon_standard_analysis_cdwpttchi2_3m_2s_omilrec_jvertex(
        new cosmo_shape_muon_standard_analysis(
            "cosmo_shape_muon_standard_analysis_cdwpttchi2_3m_2s_omilrec_jvertex", 
            filepath, 
            suffix, 
            "CdWpTtChi2", 
            timestamp{0, 5000000}, 
            timestamp{2, 0}, 
            timestamp{-2, 0}, 
            timestamp{0, -5000000}, 
            3000.0
        )
    );
    if (!registry.book(cosmo_shape_muon_standard_analysis_cdwpttchi2_3m_2s_omilrec_jvertex)) return 1;

    std::shared_ptr<analysis_base> cosmo_shape_muon_with_neutron_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex(
        new cosmo_shape_muon_with_neutron_analysis(
            "cosmo_shape_muon_with_neutron_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex", 
            filepath, 
            suffix, 
            "CdWpTtChi2", 
            timestamp{0, 5000000}, 
            timestamp{0, 1200000000}, 
            timestamp{0, -1200000000}, 
            timestamp{0, -5000000}, 
            3000.0
        )
    );
    if (!registry.book(cosmo_shape_muon_with_neutron_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex)) return 1;

    std::shared_ptr<analysis_base> cosmo_shape_muon_with_neutron_analysis_cdwpttchi2_1_5m_1_2s_omilrec_jvertex(
        new cosmo_shape_muon_with_neutron_analysis(
            "cosmo_shape_muon_with_neutron_analysis_cdwpttchi2_1_5m_1_2s_omilrec_jvertex", 
            filepath, 
            suffix, 
            "CdWpTtChi2", 
            timestamp{0, 5000000}, 
            timestamp{0, 1200000000}, 
            timestamp{0, -1200000000}, 
            timestamp{0, -5000000}, 
            1500.0
        )
    );
    if (!registry.book(cosmo_shape_muon_with_neutron_analysis_cdwpttchi2_1_5m_1_2s_omilrec_jvertex)) return 1;

    std::shared_ptr<analysis_base> cosmo_shape_muon_changing_veto_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex(
        new cosmo_shape_muon_changing_veto_analysis(
            "cosmo_shape_muon_changing_veto_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex", 
            filepath, 
            suffix, 
            "CdWpTtChi2", 
            timestamp{0, 5000000}, 
            timestamp{0, 1200000000}, 
            timestamp{0, -1200000000}, 
            timestamp{0, -5000000}, 
            3000.0
        )
    );
    if (!registry.book(cosmo_shape_muon_changing_veto_analysis_cdwpttchi2_3m_1_2s_omilrec_jvertex)) return 1;



    std::shared_ptr<analysis_base> cosmo_shape_neutron_analysis_3m_1_2s_omilred_jvertex(
        new cosmo_shape_neutron_analysis(
            "cosmo_shape_neutron_analysis_3m_1_2s_omilred_jvertex", 
            filepath, 
            suffix, 
            timestamp{0, 5000000}, 
            timestamp{0, 1200000000}, 
            timestamp{0, -1200000000}, 
            timestamp{0, -5000000}, 
            4000.0
        )
    );
    if (!registry.book(cosmo_shape_neutron_analysis_3m_1_2s_omilred_jvertex)) return 1;


    
    if (!manager.run()) return 1;
    if (!manager.result()) return 1;
    if (!manager.save()) return 1;

    return 0;
}