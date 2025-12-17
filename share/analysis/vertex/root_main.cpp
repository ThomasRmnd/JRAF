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

#include "analysis.hpp"

void daq_time(const std::string& filename) {
    TChain* chain = new TChain("DAQTree");
    if (!chain) {
        std::cerr << "Cannot create TChain DAQTree\n";
        return;
    }
    chain->Add(filename.c_str());
    time_t daq_sec;
    int daq_nsec;
    time_t muveto_sec;
    int muveto_nsec;
    chain->SetBranchAddress("daq_sec", &daq_sec);
    chain->SetBranchAddress("daq_nsec", &daq_nsec);
    chain->SetBranchAddress("muveto_sec", &muveto_sec);
    chain->SetBranchAddress("muveto_nsec", &muveto_nsec);
    timestamp tot_ts, tot_ts_mu;
    for (int k = 0; k < chain->GetEntries(); ++k) {
        chain->GetEntry(k);
        timestamp ts{daq_sec, daq_nsec};
        timestamp ts_mu{muveto_sec, muveto_nsec};
        tot_ts += ts;
        tot_ts_mu += ts_mu;
    }
    std::cout << "Total DAQ time: " << tot_ts << '\n';
    std::cout << "Total MuVeto time: " << tot_ts_mu << '\n';
    double tot_seconds = timestamp_to_double(tot_ts);
    double tot_seconds_mu = timestamp_to_double(tot_ts_mu);
    std::cout << "Total DAQ time in days: " << tot_seconds / (3600.0 * 24.0) << '\n';
    std::cout << "Total MuVeto time in days: " << tot_seconds_mu / (3600.0 * 24.0) << '\n';
}

int root_main(const std::string& filepath, const std::string& suffix) {
    daq_time(filepath);

    analysis_registry registry;
    analysis_manager manager(registry);

    std::shared_ptr<analysis_base> main_analysis(new ibd_analysis(filepath, suffix));
    if (!registry.book(main_analysis)) return 1;


    std::shared_ptr<analysis_base> cosmo_rate_with_neutron_analysis(new cosmo_rate_analysis(filepath, suffix));
    if (!registry.book(cosmo_rate_with_neutron_analysis)) return 1;

    std::shared_ptr<analysis_base> cosmo_shape_analysis_before_after(new cosmo_shape_analysis(filepath, suffix, timestamp{0, 5000000}, timestamp{0, 1200000000}, timestamp{0, -1200000000}, timestamp{0, -5000000}, 3000.0));
    std::shared_ptr<analysis_base> cosmo_shape_analysis_after_later(new cosmo_shape_analysis(filepath, suffix, timestamp{0, 5000000}, timestamp{0, 1200000000}, timestamp{0, 5005000000}, timestamp{0, 6200000000}, 3000.0));
    if (!registry.book(cosmo_shape_analysis_before_after)) return 1;
    if (!registry.book(cosmo_shape_analysis_after_later)) return 1;
    
    if (!manager.run()) return 1;

    return 0;
}