#include <iostream>

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
    h->GetXaxis()->Center(true);
    h->GetYaxis()->Center(true);
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
}

int fast_muon_rate_estimation(const char* filepath) {
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
    std::vector<std::string> *method;
    std::vector<unsigned char> *det;
    std::vector<double> *quality;
    std::vector<double> *iposx, *iposy, *iposz, *fposx, *fposy, *fposz;
    tree->SetBranchAddress("run_id", &run_id);
    tree->SetBranchAddress("sec", &sec);
    tree->SetBranchAddress("nsec", &nsec);
    tree->SetBranchAddress("totq_cd", &totq_cd);
    tree->SetBranchAddress("totq_wp", &totq_wp);
    tree->SetBranchAddress("method", &method);
    tree->SetBranchAddress("det", &det);
    tree->SetBranchAddress("quality", &quality);
    tree->SetBranchAddress("iposx", &iposx);
    tree->SetBranchAddress("iposy", &iposy);
    tree->SetBranchAddress("iposz", &iposz);
    tree->SetBranchAddress("fposx", &fposx);
    tree->SetBranchAddress("fposy", &fposy);
    tree->SetBranchAddress("fposz", &fposz);

    TH1D* h_time_to_previous_muon_cd = new TH1D("h_time_to_previous_muon_cd", "Time to previous muon in CD; #Delta t (s); Entries;", 100, 0.0, 2.0);
    TH1D* h_time_to_previous_muon_wp = new TH1D("h_time_to_previous_muon_wp", "Time to previous muon in WP; #Delta t (s); Entries;", 100, 0.0, 2.0);
    TH1D* h_time_to_previous_muon = new TH1D("h_time_to_previous_muon", "Time to previous muon; #Delta t (s); Entries;", 100, 0.0, 2.0);

    TTimeStamp prvts_cd{0, 0}, prvts_wp{0, 0}, prvts{0, 0};
    long nentries = tree->GetEntries();
    for (long k = 0l; k < nentries; ++k) {
        tree->GetEntry(k);
        TTimeStamp ts{sec, nsec};
        if (totq_cd > 0.0 && totq_wp > 0.0) {
            h_time_to_previous_muon_cd->Fill(ts - prvts_cd);
            h_time_to_previous_muon_wp->Fill(ts - prvts_wp);
            h_time_to_previous_muon->Fill(ts - prvts);
            prvts_cd = ts;
            prvts_wp = ts;
            prvts = ts;
        }
        else if (totq_cd <= 0.0 && totq_wp > 0.0) {
            h_time_to_previous_muon_wp->Fill(ts - prvts_wp);
            prvts_wp = ts;
        }
        else if (totq_cd > 0.0 && totq_wp <= 0.0) {
            h_time_to_previous_muon_cd->Fill(ts - prvts_cd);
            prvts_cd = ts;
        }
    }

    if (int res = fit_and_plot_rate(h_time_to_previous_muon_cd)) return res;
    if (int res = fit_and_plot_rate(h_time_to_previous_muon_wp)) return res;
    if (int res = fit_and_plot_rate(h_time_to_previous_muon)) return res;

    return 0;
}