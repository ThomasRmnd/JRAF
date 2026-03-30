#include <cstdint>
#include <iostream>
#include <limits>
#include <unordered_map>

#include <TFile.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TVector3.h>

TF1* fit_rate(TH1D* h) {
    TF1* f = new TF1(Form("f_%s", h->GetName()), "[0] * exp(-[1] * x)", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    f->SetParameter(0, h->GetBinContent(1));
    f->SetParameter(1, h->GetRMS());
    h->Fit(f, "R");
    return f;
}

int fast_muon_classification_comparison(const char* filepath) {
    TChain* chain = new TChain("muons");
    if (!chain) {
        std::cerr << "Cannot create chain muons\n";
        return 1;
    }
    chain->Add(filepath);

    int run_id;
    time_t sec;
    int nsec;
    double totq_cd;
    double totq_wp;
    std::vector<std::string> *method = nullptr;
    std::vector<unsigned char> *det = nullptr;
    std::vector<double> *quality = nullptr;
    std::vector<double> *iposx = nullptr, *iposy = nullptr, *iposz = nullptr; 
    std::vector<double> *fposx = nullptr, *fposy = nullptr, *fposz = nullptr;
    chain->SetBranchAddress("run_id", &run_id);
    chain->SetBranchAddress("sec", &sec);
    chain->SetBranchAddress("nsec", &nsec);
    chain->SetBranchAddress("totq_cd", &totq_cd);
    chain->SetBranchAddress("totq_wp", &totq_wp);
    chain->SetBranchAddress("method", &method);
    chain->SetBranchAddress("det", &det);
    chain->SetBranchAddress("quality", &quality);
    chain->SetBranchAddress("iposx", &iposx);
    chain->SetBranchAddress("iposy", &iposy);
    chain->SetBranchAddress("iposz", &iposz);
    chain->SetBranchAddress("fposx", &fposx);
    chain->SetBranchAddress("fposy", &fposy);
    chain->SetBranchAddress("fposz", &fposz);

    std::unordered_map<int, TH1D*> h_time_to_previous_muon_single;
    std::unordered_map<int, TH1D*> h_time_to_previous_muon_bundle;
    std::unordered_map<int, TTimeStamp> prvts_single;
    std::unordered_map<int, TTimeStamp> prvts_bundle;
    std::unordered_map<int, TF1*> fit_res_single;
    std::unordered_map<int, TF1*> fit_res_bundle;

    int min_run_id = std::numeric_limits<int>::max();
    int max_run_id = std::numeric_limits<int>::min();

    Long64_t nentries = chain->GetEntries();
    std::cout << "Info: Found " << nentries << " entries in joint reco files\n";
    for (Long64_t k = 0l; k < nentries; ++k) {
        chain->GetEntry(k);
        min_run_id = std::min(min_run_id, run_id);
        max_run_id = std::max(max_run_id, run_id);
        std::size_t ntracks_cdclassify = 0ul;
        std::size_t ntracks_wpclassify = 0ul;
        for (std::size_t i = 0ul; i < method->size(); ++i) {
            if ((*method)[i] == "CdClassify") {
                ++ntracks_cdclassify;
            }
            if ((*method)[i] == "WpBasic") {
                ++ntracks_wpclassify;
            }
        }
        if (ntracks_wpclassify == 1) {
            if (prvts_single.find(run_id) == prvts_single.end()) {
                prvts_single[run_id] = TTimeStamp{sec, nsec};
                h_time_to_previous_muon_single[run_id] = new TH1D(Form("h_time_to_previous_muon_single_%d", run_id), Form("Time to previous muon for run %d; #Delta t (s); Entries;", run_id), 100, 0.0, 5.0);
                continue;
            }
            TTimeStamp ts{sec, nsec};
            h_time_to_previous_muon_single[run_id]->Fill(ts - prvts_single[run_id]);
            prvts_single[run_id] = ts;
        }
        if (ntracks_wpclassify > 1) {
            if (prvts_bundle.find(run_id) == prvts_bundle.end()) {
                prvts_bundle[run_id] = TTimeStamp{sec, nsec};
                h_time_to_previous_muon_bundle[run_id] = new TH1D(Form("h_time_to_previous_muon_bundle_%d", run_id), Form("Time to previous muon for run %d; #Delta t (s); Entries;", run_id), 100, 0.0, 5.0);
                continue;
            }
            TTimeStamp ts{sec, nsec};
            h_time_to_previous_muon_bundle[run_id]->Fill(ts - prvts_bundle[run_id]);
            prvts_bundle[run_id] = ts;
        }
    }

    for (const auto& [run_id, h] : h_time_to_previous_muon_single) {
        fit_res_single[run_id] = fit_rate(h);
    }
    for (const auto& [run_id, h] : h_time_to_previous_muon_bundle) {
        fit_res_bundle[run_id] = fit_rate(h);
    }

    TH1D* h_rate_per_run_single = new TH1D("h_rate_per_run_single", "Rate per run for single;RUN ID;Rate (cps);", max_run_id - min_run_id + 1, min_run_id, max_run_id + 1);
    TH1D* h_rate_per_run_bundle = new TH1D("h_rate_per_run_bundle", "Rate per run for bundle;RUN ID;Rate (cps);", max_run_id - min_run_id + 1, min_run_id, max_run_id + 1);
    for (const auto& [run_id, h] : h_time_to_previous_muon_single) {
        h_rate_per_run_single->SetBinContent(run_id - min_run_id + 1, fit_res_single[run_id]->GetParameter(1));
        h_rate_per_run_single->SetBinError(run_id - min_run_id + 1, fit_res_single[run_id]->GetParError(1));
    }
    for (const auto& [run_id, h] : h_time_to_previous_muon_bundle) {
        h_rate_per_run_bundle->SetBinContent(run_id - min_run_id + 1, fit_res_bundle[run_id]->GetParameter(1));
        h_rate_per_run_bundle->SetBinError(run_id - min_run_id + 1, fit_res_bundle[run_id]->GetParError(1));
    }

    TCanvas* c_rate_per_run = new TCanvas("c_rate_per_run", "Rate per run", 1000, 1000);
    c_rate_per_run->cd();
    h_rate_per_run_single->SetStats(0);
    h_rate_per_run_single->GetYaxis()->SetMaxDigits(3);
    h_rate_per_run_single->GetXaxis()->CenterTitle(true);
    h_rate_per_run_single->GetYaxis()->CenterTitle(true);
    h_rate_per_run_single->GetYaxis()->SetTitleOffset(1.25);
    h_rate_per_run_single->SetMarkerStyle(kFullCircle);
    h_rate_per_run_single->SetMarkerColor(kBlue);
    h_rate_per_run_single->SetMarkerSize(1.0);
    h_rate_per_run_single->SetLineWidth(2);
    h_rate_per_run_single->SetLineColor(kBlue);
    h_rate_per_run_single->SetMinimum(0.0);
    h_rate_per_run_single->Draw();
    c_rate_per_run->Update();
    h_rate_per_run_bundle->SetStats(0);
    h_rate_per_run_bundle->GetYaxis()->SetMaxDigits(3);
    h_rate_per_run_bundle->GetXaxis()->CenterTitle(true);
    h_rate_per_run_bundle->GetYaxis()->CenterTitle(true);
    h_rate_per_run_bundle->GetYaxis()->SetTitleOffset(1.25);
    h_rate_per_run_bundle->SetMarkerStyle(kFullCircle);
    h_rate_per_run_bundle->SetMarkerColor(kRed);
    h_rate_per_run_bundle->SetMarkerSize(1.0);
    h_rate_per_run_bundle->SetLineWidth(2);
    h_rate_per_run_bundle->SetLineColor(kRed);
    h_rate_per_run_bundle->SetMinimum(0.0);
    h_rate_per_run_bundle->Draw("SAME");
    c_rate_per_run->Update();
    c_rate_per_run->SetTickx();
    c_rate_per_run->SetTicky();
    c_rate_per_run->SetGrid();
    c_rate_per_run->Update();

    return 0;
}