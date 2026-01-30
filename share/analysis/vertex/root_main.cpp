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

struct veto_info {
    int run_id;
    unsigned char veto_type;
    time_t veto_sec;
    int veto_nsec;
};

void daq(const std::string& filename) {
    TChain* chain_daq = new TChain("DAQ");
    TChain* chain_veto = new TChain("Veto");
    if (!chain_daq || !chain_veto) {
        std::cerr << "Cannot create TChain DAQ or Veto\n";
        return;
    }
    
    chain_daq->Add(filename.c_str());
    DAQ daq;
    chain_daq->SetBranchAddress("run_id", &daq.run_id);
    chain_daq->SetBranchAddress("sec", &daq.sec);
    chain_daq->SetBranchAddress("nsec", &daq.nsec);
    std::map<int, timestamp> daq_map;
    timestamp tot_ts;

    chain_veto->Add(filename.c_str());
    Veto veto;
    chain_veto->SetBranchAddress("run_id", &veto.run_id);
    chain_veto->SetBranchAddress("sec", &veto.sec);
    chain_veto->SetBranchAddress("nsec", &veto.nsec);
    chain_veto->SetBranchAddress("veto_type", &veto.veto_type);
    chain_veto->SetBranchAddress("veto_sec", &veto.veto_sec);
    chain_veto->SetBranchAddress("veto_nsec", &veto.veto_nsec);
    std::map<timestamp, veto_info> veto_map;
    std::map<int, std::map<unsigned char, std::size_t>> veto_type_per_run_map;

    for (long k = 0; k < chain_daq->GetEntries(); ++k) {
        chain_daq->GetEntry(k);
        timestamp ts{daq.sec, daq.nsec};
        daq_map[daq.run_id] += ts;
        tot_ts += ts;
    }

    for (long k = 0; k < chain_veto->GetEntries(); ++k) {
        chain_veto->GetEntry(k);
        timestamp ts{veto.sec, veto.nsec};
        veto_map[ts] = {veto.run_id, veto.veto_type, veto.veto_sec, veto.veto_nsec};
    }

    timestamp prev_mu_ts;
    int min_run_id = daq_map.size() ? daq_map.begin()->first : 0;
    int max_run_id = daq_map.size() ? (--daq_map.end())->first : 0;
    TH1D* h_mu_tot_rate = new TH1D("h_mu_tot_rate", "h_mu_tot_rate", 100, 0.0, 5.0);
    for (std::map<timestamp, veto_info>::iterator it = veto_map.begin(); it != veto_map.end(); ++it) {
        veto_type_per_run_map[it->second.run_id][it->second.veto_type] += 1;
        if (prev_mu_ts == timestamp{0, 0}) {
            prev_mu_ts = it->first;
            continue;
        }
        h_mu_tot_rate->Fill(timestamp_to_double(it->first - prev_mu_ts));
        prev_mu_ts = it->first;
    }
    
    std::cout << "Total DAQ time: " << tot_ts << '\n';
    double tot_seconds = timestamp_to_double(tot_ts);
    std::cout << "Total DAQ time in days: " << tot_seconds / (3600.0 * 24.0) << '\n';

    TCanvas* c_mu_tot_rate = new TCanvas("c_mu_tot_rate", "c_mu_tot_rate", 1000, 1000);
    c_mu_tot_rate->cd();
    h_mu_tot_rate->Draw("HIST");
    c_mu_tot_rate->SetLogy();
    c_mu_tot_rate->SetGrid();
    c_mu_tot_rate->SetTickx();
    c_mu_tot_rate->SetTicky();
    c_mu_tot_rate->Update();

    TH1D* h_beg_veto_per_run = new TH1D("h_beg_veto_per_run", "h_beg_veto_per_run", max_run_id - min_run_id + 1, min_run_id, max_run_id);
    TH1D* h_hdr_veto_per_run = new TH1D("h_hdr_veto_per_run", "h_hdr_veto_per_run", max_run_id - min_run_id + 1, min_run_id, max_run_id);
    TH1D* h_gap_veto_per_run = new TH1D("h_gap_veto_per_run", "h_gap_veto_per_run", max_run_id - min_run_id + 1, min_run_id, max_run_id);
    TH1D* h_mu_veto_per_run = new TH1D("h_mu_veto_per_run", "h_mu_veto_per_run", max_run_id - min_run_id + 1, min_run_id, max_run_id);
    for (std::map<int, std::map<unsigned char, std::size_t>>::iterator it = veto_type_per_run_map.begin(); it != veto_type_per_run_map.end(); ++it) {
        for (std::map<unsigned char, std::size_t>::iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
            if (jt->first == 1) {
                h_beg_veto_per_run->SetBinContent(it->first - min_run_id + 1, jt->second);
            }
            else if (jt->first == 2) {
                h_hdr_veto_per_run->SetBinContent(it->first - min_run_id + 1, jt->second);
            }
            else if (jt->first == 3) {
                h_gap_veto_per_run->SetBinContent(it->first - min_run_id + 1, jt->second);
            }
            else if (jt->first == 4) {
                h_mu_veto_per_run->SetBinContent(it->first - min_run_id + 1, jt->second);
            }
        }
    }

    TCanvas* c_veto_per_run = new TCanvas("c_veto_per_run", "c_veto_per_run", 1000, 1000);
    c_veto_per_run->cd();
    h_beg_veto_per_run->SetStats(false);
    h_beg_veto_per_run->SetMaximum(std::max({h_beg_veto_per_run->GetMaximum(), h_hdr_veto_per_run->GetMaximum(), h_gap_veto_per_run->GetMaximum(), h_mu_veto_per_run->GetMaximum()}));
    h_beg_veto_per_run->GetXaxis()->SetTitle("Run ID");
    h_beg_veto_per_run->GetYaxis()->SetTitle("Entries");
    h_beg_veto_per_run->Draw("HIST");
    h_hdr_veto_per_run->Draw("SAME");
    h_gap_veto_per_run->Draw("SAME");
    h_mu_veto_per_run->Draw("SAME");
    c_veto_per_run->SetLogy();
    c_veto_per_run->SetGrid();
    c_veto_per_run->SetTickx();
    c_veto_per_run->SetTicky();
    c_veto_per_run->Update();
}

int root_main(const std::string& filepath) {
    std::string suffix = "__OMILREC_JVtx";

    daq(filepath);

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