#include <iostream>
#include <limits>
#include <unordered_map>

#include <TFile.h>
#include <TTimeStamp.h>
#include <TTree.h>

int fit_and_plot_rate(TH1D* h) {
    TF1* f = new TF1(Form("f_%s", h->GetName()), "[0] * exp(-[1] * x)", 0.0, 2.0);
    
    TCanvas* c = new TCanvas(Form("c_%s", h->GetName()), "Time to previous muon", 1000, 1000);
    c->cd();

    f->SetParameter(0, h->GetBinContent(1));
    f->SetParameter(1, h->GetRMS());
    h->Fit(f, "R");
    
    h->SetLineColor(kBlue);
    h->SetLineWidth(3);
    h->GetXaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetMaxDigits(3);
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->CenterTitle(true);
    h->GetYaxis()->SetTitleOffset(1.25);
    h->Draw();
    f->SetLineColor(kRed);
    f->SetLineWidth(3);
    f->Draw("SAME");
    
    c->Update();

    TPaveStats* ps = (TPaveStats*)h->FindObject("stats");
    ps->SetOptStat(1110);
    ps->SetOptFit(1112); 
    ps->SetX1NDC(0.50);
    ps->SetX2NDC(0.85);
    ps->SetY1NDC(0.50);
    ps->SetY2NDC(0.85);

    c->Modified();
    c->SetTickx();
    c->SetTicky();
    c->SetGrid();
    c->Update();

    return 0;
}

TF1* fit_rate(TH1D* h) {
    TF1* f = new TF1(Form("f_%s", h->GetName()), "[0] * exp(-[1] * x)", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    f->SetParameter(0, h->GetBinContent(1));
    f->SetParameter(1, h->GetRMS());
    h->Fit(f, "R");
    return f;
}

int fast_muon_rate_estimation_cdwpttchi2(const char* filepath) {
    TFile* file = TFile::Open(filepath, "READ");
    if (!file) {
        std::cerr << "Cannot open file " << filepath << '\n';
        return 1;
    }
    TTree* tree = file->Get<TTree>("muons");
    if (!tree) {
        std::cerr << "Cannot retrieve tree muons in file " << filepath << '\n';
        return 1;
    }

    int run_id;
    time_t sec;
    int nsec;
    double totq_cd;
    double totq_wp;
    double chi2;
    double iposx, iposy, iposz; 
    double fposx, fposy, fposz;
    tree->SetBranchAddress("run_id", &run_id);
    tree->SetBranchAddress("sec", &sec);
    tree->SetBranchAddress("nsec", &nsec);
    tree->SetBranchAddress("totq_cd", &totq_cd);
    tree->SetBranchAddress("totq_wp", &totq_wp);
    tree->SetBranchAddress("chi2", &chi2);
    tree->SetBranchAddress("iposx", &iposx);
    tree->SetBranchAddress("iposy", &iposy);
    tree->SetBranchAddress("iposz", &iposz);
    tree->SetBranchAddress("fposx", &fposx);
    tree->SetBranchAddress("fposy", &fposy);
    tree->SetBranchAddress("fposz", &fposz);

    std::unordered_map<int, TH1D*> h_time_to_previous_muon;
    std::unordered_map<int, TTimeStamp> prvts;
    std::unordered_map<int, TF1*> fit_res;

    // TH1D* h_time_to_previous_muon_cdwpttchi2 = new TH1D("h_time_to_previous_muon_cdwpttchi2", "Time to previous muon for CdWpTtChi2; #Delta t (s); Entries;", 100, 0.0, 5.0);
    // TTimeStamp prvts{0, 0};

    int min_run_id = std::numeric_limits<int>::max();
    int max_run_id = std::numeric_limits<int>::min();

    Long64_t nentries = tree->GetEntries();
    for (Long64_t k = 0l; k < nentries; ++k) {
        tree->GetEntry(k);
        min_run_id = std::min(min_run_id, run_id);
        max_run_id = std::max(max_run_id, run_id);
        if (prvts.find(run_id) == prvts.end()) {
            prvts[run_id] = TTimeStamp{sec, nsec};
            h_time_to_previous_muon[run_id] = new TH1D(Form("h_time_to_previous_muon_%d", run_id), Form("Time to previous muon for run %d; #Delta t (s); Entries;", run_id), 100, 0.0, 5.0);
            continue;
        }
        TTimeStamp ts{sec, nsec};
        h_time_to_previous_muon[run_id]->Fill(ts - prvts[run_id]);
        prvts[run_id] = ts;
    }

    for (const auto& [run_id, h] : h_time_to_previous_muon) {
        fit_res[run_id] = fit_rate(h);
    }

    TH1D* h_rate_per_run = new TH1D("h_rate_per_run", "Rate per run;RUN ID;Rate (cps);", max_run_id - min_run_id + 1, min_run_id, max_run_id + 1);
    for (const auto& [run_id, h] : h_time_to_previous_muon) {
        h_rate_per_run->SetBinContent(run_id - min_run_id + 1, fit_res[run_id]->GetParameter(1));
        h_rate_per_run->SetBinError(run_id - min_run_id + 1, fit_res[run_id]->GetParError(1));
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

    // if (int res = fit_and_plot_rate(h_time_to_previous_muon_cdwpttchi2)) return res;

    return 0;
}