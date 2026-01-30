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
#include "analysis/cosmo_shape_analysis.hpp"
#include "analysis/cosmo_shape_neutron_analysis.hpp"
#include "analysis/ibd_analysis.hpp"
#include "analysis/ibd_muon_veto_analysis.hpp"

struct DAQ {
    int run_id;
    time_t sec;
    int nsec;
};

struct Veto {
    int run_id;
    time_t sec;
    int nsec;
    unsigned char veto_type;
    time_t veto_sec;
    int veto_nsec;
};

struct MuonInfo {
    int run_id;
    time_t sec;
    int nsec;
    float totq_cd;
    float totq_wp;
    unsigned char det;
};

void save_meta_info(const std::string& filename) {
    TChain* chain_daq = new TChain("DAQ");
    TChain* chain_veto = new TChain("Veto");
    TChain* chain_muon = new TChain("MuonInfo");
    if (!chain_daq || !chain_veto || !chain_muon) {
        std::cerr << "Cannot create TChain DAQ or Veto or MuonInfo\n";
        return;
    }
    
    chain_daq->Add(filename.c_str());
    DAQ daq;
    chain_daq->SetBranchAddress("run_id", &daq.run_id);
    chain_daq->SetBranchAddress("sec", &daq.sec);
    chain_daq->SetBranchAddress("nsec", &daq.nsec);

    chain_veto->Add(filename.c_str());
    Veto veto;
    chain_veto->SetBranchAddress("run_id", &veto.run_id);
    chain_veto->SetBranchAddress("sec", &veto.sec);
    chain_veto->SetBranchAddress("nsec", &veto.nsec);
    chain_veto->SetBranchAddress("veto_type", &veto.veto_type);
    chain_veto->SetBranchAddress("veto_sec", &veto.veto_sec);
    chain_veto->SetBranchAddress("veto_nsec", &veto.veto_nsec);

    chain_muon->Add(filename.c_str());
    MuonInfo muon;
    chain_muon->SetBranchAddress("run_id", &muon.run_id);
    chain_muon->SetBranchAddress("sec", &muon.sec);
    chain_muon->SetBranchAddress("nsec", &muon.nsec);
    chain_muon->SetBranchAddress("totq_cd", &muon.totq_cd);
    chain_muon->SetBranchAddress("totq_wp", &muon.totq_wp);
    chain_muon->SetBranchAddress("det", &muon.det);

    TFile* f_run_info = TFile::Open("run_info.root", "RECREATE");
    if (!f_run_info) {
        std::cerr << "Cannot open file run_info.root for writing\n";
        return;
    }
    f_run_info->cd();

    TTree* out_daq = chain_daq->CloneTree(0);
    for (Long64_t i = 0; i < chain_daq->GetEntries(); ++i) {
        chain_daq->GetEntry(i);
        out_daq->Fill();
    }
    out_daq->Write();

    TTree* out_veto = chain_veto->CloneTree(0);
    for (Long64_t i = 0; i < chain_veto->GetEntries(); ++i) {
        chain_veto->GetEntry(i);
        out_veto->Fill();
    }
    out_veto->Write();

    TTree* out_muon = chain_muon->CloneTree(0);
    for (Long64_t i = 0; i < chain_muon->GetEntries(); ++i) {
        chain_muon->GetEntry(i);
        out_muon->Fill();
    }
    out_muon->Write();

    f_run_info->Write(); 
    f_run_info->Close();
    
    std::cout << "Successfully saved meta info to run_info.root" << std::endl;
}

int root_main(const std::string& filepath) {
    std::string suffix = "__OMILREC_JVtx";

    save_meta_info(filepath);

    analysis_registry registry;
    analysis_manager manager(registry);

    std::shared_ptr<analysis_base> main_analysis(new ibd_analysis("ibd_analysis", filepath, suffix));
    if (!registry.book(main_analysis)) return 1;

    std::shared_ptr<analysis_base> ibd_with_muon_veto_analysis(new ibd_muon_veto_analysis("ibd_muon_veto_analysis_cdwpttchi2", filepath, suffix, "CdWpTtChi2", timestamp{0, 5000000}, timestamp{0, 1200000000}, 3000.0));
    if (!registry.book(ibd_with_muon_veto_analysis)) return 1;

    std::shared_ptr<analysis_base> cosmo_rate_with_neutron_analysis(new cosmo_rate_analysis("cosmo_rate_analysis", filepath, suffix));
    // if (!registry.book(cosmo_rate_with_neutron_analysis)) return 1;

    std::shared_ptr<analysis_base> cosmo_shape_analysis_before_after_cdwpttchi2(new cosmo_shape_analysis("cosmo_shape_analysis_cdwpttchi2", filepath, suffix, "CdWpTtChi2", timestamp{0, 5000000}, timestamp{0, 1200000000}, timestamp{0, -1200000000}, timestamp{0, -5000000}, 3000.0));
    if (!registry.book(cosmo_shape_analysis_before_after_cdwpttchi2)) return 1;

    std::shared_ptr<analysis_base> cosmo_shape_analysis_before_after_cdclassify(new cosmo_shape_analysis("cosmo_shape_analysis_cdclassify", filepath, suffix, "CdClassify", timestamp{0, 5000000}, timestamp{0, 1200000000}, timestamp{0, -1200000000}, timestamp{0, -5000000}, 3000.0));
    // if (!registry.book(cosmo_shape_analysis_before_after_cdclassify)) return 1;

    std::shared_ptr<analysis_base> cosmo_shape_analysis_before_after_tt(new cosmo_shape_analysis("cosmo_shape_analysis_tt", filepath, suffix, "Tt", timestamp{0, 5000000}, timestamp{0, 1200000000}, timestamp{0, -1200000000}, timestamp{0, -5000000}, 3000.0));
    // if (!registry.book(cosmo_shape_analysis_before_after_tt)) return 1;

    std::shared_ptr<analysis_base> cosmo_shape_with_neutron_analysis(new cosmo_shape_neutron_analysis("cosmo_shape_neutron_analysis", filepath, suffix, timestamp{0, 5000000}, timestamp{0, 1200000000}, timestamp{0, -1200000000}, timestamp{0, -5000000}, 4000.0));
    if (!registry.book(cosmo_shape_with_neutron_analysis)) return 1;
    
    if (!manager.run()) return 1;
    // if (!manager.result()) return 1;
    if (!manager.save()) return 1;

    return 0;
}