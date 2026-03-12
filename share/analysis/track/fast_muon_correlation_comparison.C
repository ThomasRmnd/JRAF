#include <cstdint>
#include <iostream>
#include <limits>
#include <set>
#include <unordered_map>

#include <TFile.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TVector3.h>

struct track {

    int run_id;
    TTimeStamp ts;

};

bool operator<(const track& t1, const track& t2) {
    return t1.ts < t2.ts;
}

TF1* fit_rate(TH1D* h) {
    TF1* f = new TF1(Form("f_%s", h->GetName()), "[0] * exp(-[1] * x)", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    f->SetParameter(0, h->GetBinContent(1));
    f->SetParameter(1, h->GetRMS());
    h->Fit(f, "R");
    return f;
}

int fast_muon_correlation_comparison(const char* filepath_cdwpttchi2, const char* filepath_amber) {
    TChain* chain_cdwpttchi2 = new TChain("muons");
    if (!chain_cdwpttchi2) {
        std::cerr << "Cannot create chain muons\n";
        return 1;
    }
    chain_cdwpttchi2->Add(filepath_cdwpttchi2);

    int run_id;
    time_t sec;
    int nsec;
    double totq_cd;
    double totq_wp;
    double chi2;
    double iposx, iposy, iposz; 
    double fposx, fposy, fposz;
    unsigned int ntracks_cdclassify, ntracks_wpclassify;
    unsigned int nstoppings_cdclassify, nstoppings_wpclassify;
    chain_cdwpttchi2->SetBranchAddress("run_id", &run_id);
    chain_cdwpttchi2->SetBranchAddress("sec", &sec);
    chain_cdwpttchi2->SetBranchAddress("nsec", &nsec);
    chain_cdwpttchi2->SetBranchAddress("totq_cd", &totq_cd);
    chain_cdwpttchi2->SetBranchAddress("totq_wp", &totq_wp);
    chain_cdwpttchi2->SetBranchAddress("chi2", &chi2);
    chain_cdwpttchi2->SetBranchAddress("iposx", &iposx);
    chain_cdwpttchi2->SetBranchAddress("iposy", &iposy);
    chain_cdwpttchi2->SetBranchAddress("iposz", &iposz);
    chain_cdwpttchi2->SetBranchAddress("fposx", &fposx);
    chain_cdwpttchi2->SetBranchAddress("fposy", &fposy);
    chain_cdwpttchi2->SetBranchAddress("fposz", &fposz);
    chain_cdwpttchi2->SetBranchAddress("ntracks_cdclassify", &ntracks_cdclassify);
    chain_cdwpttchi2->SetBranchAddress("ntracks_wpclassify", &ntracks_wpclassify);
    chain_cdwpttchi2->SetBranchAddress("nstoppings_cdclassify", &nstoppings_cdclassify);
    chain_cdwpttchi2->SetBranchAddress("nstoppings_wpclassify", &nstoppings_wpclassify);

    std::<set> tracks_cdwpttchi2;
    Long64_t nentries = chain_cdwpttchi2->GetEntries();
    for (Long64_t k = 0l; k < nentries; ++k) {
        chain_cdwpttchi2->GetEntry(k);
        tracks_cdwpttchi2.insert(track{run_id, TTimeStamp{sec, nsec}});
    }

    TChain* chain_amber = new TChain("MuonReco");
    if (!chain_amber) {
        std::cerr << "Cannot create chain MuonReco\n";
        return 1;
    }
    chain_amber->Add(filepath_amber);

    int runID;
    int eventID;
    int fSec;
    int fNanoSec;
    int muonType;
    float xin, yin, zin;
    float xout, yout, zout;
    float charge;
    chain->SetBranchAddress("runID", &runID);
    chain->SetBranchAddress("eventID", &eventID);
    chain->SetBranchAddress("fSec", &fSec);
    chain->SetBranchAddress("fNanoSec", &fNanoSec);
    chain->SetBranchAddress("muonType", &muonType);
    chain->SetBranchAddress("xin", &xin);
    chain->SetBranchAddress("yin", &yin);
    chain->SetBranchAddress("zin", &zin);
    chain->SetBranchAddress("xout", &xout);
    chain->SetBranchAddress("yout", &yout);
    chain->SetBranchAddress("zout", &zout);
    chain->SetBranchAddress("charge", &charge);

    std::set<track> tracks_amber;
    nentries = chain_amber->GetEntries();
    for (Long64_t k = 0l; k < nentries; ++k) {
        chain_amber->GetEntry(k);
        tracks_amber.insert(track{runID, TTimeStamp{fSec, fNanoSec}});
    }

    // Find a +/- 1000 ns correlation for each muons

    std::unordered_map<int, TH1D*> h_time_to_previous_muon;
    std::unordered_map<int, TTimeStamp> prvts;
    std::unordered_map<int, TF1*> fit_res;

    int min_run_id = std::numeric_limits<int>::max();
    int max_run_id = std::numeric_limits<int>::min();

    for (std::set<track>::const_iterator it = tracks_cdwpttchi2.begin(); it != tracks_cdwpttchi2.end(); ++it) {
        TTimeStamp lower(it->ts.GetSec(), it->ts.GetNanoSec() - 500);
        TTimeStamp upper(it->ts.GetSec(), it->ts.GetNanoSec() + 500);
        std::set<track>::const_iterator jt = tracks_amber.lower_bound(track{it->run_id, lower});
        if (jt != tracks_amber.end() && jt->ts < upper) {
            min_run_id = std::min(min_run_id, it->run_id);
            max_run_id = std::max(max_run_id, it->run_id);
            if (prvts.find(it->run_id) == prvts.end()) {
                prvts[it->run_id] = TTimeStamp{sec, nsec};
                h_time_to_previous_muon[it->run_id] = new TH1D(Form("h_time_to_previous_muon_%d", it->run_id), Form("Time to previous muon for run %d; #Delta t (s); Entries;", it->run_id), 100, 0.0, 5.0);
                continue;
            }
            TTimeStamp ts{sec, nsec};
            h_time_to_previous_muon[it->run_id]->Fill(ts - prvts[it->run_id]);
            prvts[it->run_id] = ts;
        }
    }

    for (const auto& [run_id, h] : h_time_to_previous_muon) {
        fit_res_single[run_id] = fit_rate(h);
    }

    TH1D* h_rate_per_run = new TH1D("h_rate_per_run", "Rate per run;RUN ID;Rate (cps);", max_run_id - min_run_id + 1, min_run_id, max_run_id + 1);
    for (const auto& [run_id, h] : h_time_to_previous_muon) {
        h_rate_per_run->SetBinContent(run_id - min_run_id + 1, fit_res_single[run_id]->GetParameter(1));
        h_rate_per_run->SetBinError(run_id - min_run_id + 1, fit_res_single[run_id]->GetParError(1));
    }

    TCanvas* c_rate_per_run = new TCanvas("c_rate_per_run", "Rate per run", 1000, 1000);
    c_rate_per_run->cd();
    h_rate_per_run->SetStats(0);
    h_rate_per_run->GetYaxis()->SetMaxDigits(3);
    h_rate_per_run->GetXaxis()->CenterTitle(true);
    h_rate_per_run->GetYaxis()->CenterTitle(true);
    h_rate_per_run->GetYaxis()->SetTitleOffset(1.25);
    h_rate_per_run->SetMarkerStyle(kFullCircle);
    h_rate_per_run->SetMarkerColor(kBlue);
    h_rate_per_run->SetMarkerSize(1.0);
    h_rate_per_run->SetLineWidth(2);
    h_rate_per_run->SetLineColor(kBlue);
    h_rate_per_run->SetMinimum(0.0);
    h_rate_per_run->Draw();
    c_rate_per_run->Update();
    c_rate_per_run->SetTickx();
    c_rate_per_run->SetTicky();
    c_rate_per_run->SetGrid();
    c_rate_per_run->Update();

    return 0;
}