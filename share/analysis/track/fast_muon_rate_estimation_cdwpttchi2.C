#include <cmath>
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
    unsigned int ntracks_cdclassify, ntracks_wpclassify;
    unsigned int nstoppings_cdclassify, nstoppings_wpclassify;
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
    tree->SetBranchAddress("ntracks_cdclassify", &ntracks_cdclassify);
    tree->SetBranchAddress("ntracks_wpclassify", &ntracks_wpclassify);
    tree->SetBranchAddress("nstoppings_cdclassify", &nstoppings_cdclassify);
    tree->SetBranchAddress("nstoppings_wpclassify", &nstoppings_wpclassify);

    std::unordered_map<int, TH1D*> h_time_to_previous_muon;
    std::unordered_map<int, TH1D*> h_time_to_previous_muon_single;
    std::unordered_map<int, TH1D*> h_time_to_previous_muon_bundle;
    std::unordered_map<int, TTimeStamp> prvts;
    std::unordered_map<int, TTimeStamp> prvts_single;
    std::unordered_map<int, TTimeStamp> prvts_bundle;
    std::unordered_map<int, TF1*> fit_res;
    std::unordered_map<int, TF1*> fit_res_single;
    std::unordered_map<int, TF1*> fit_res_bundle;

    std::unordered_map<int, double> sum_chi2;
    std::unordered_map<int, double> sum_chi2sq;
    std::unordered_map<int, std::size_t> ntracks;
    std::unordered_map<int, TH1D*> h_clippingness;

    // TH1D* h_time_to_previous_muon_cdwpttchi2 = new TH1D("h_time_to_previous_muon_cdwpttchi2", "Time to previous muon for CdWpTtChi2; #Delta t (s); Entries;", 100, 0.0, 5.0);
    // TTimeStamp prvts{0, 0};

    int min_run_id = std::numeric_limits<int>::max();
    int max_run_id = std::numeric_limits<int>::min();

    Long64_t nentries = tree->GetEntries();
    for (Long64_t k = 0l; k < nentries; ++k) {
        tree->GetEntry(k);
        min_run_id = std::min(min_run_id, run_id);
        max_run_id = std::max(max_run_id, run_id);
        if (!std::isnan(chi2) && !std::isinf(chi2)) {
            if (sum_chi2.find(run_id) == sum_chi2.end()) {
                sum_chi2[run_id] = 0.0;
                sum_chi2sq[run_id] = 0.0;
                ntracks[run_id] = 0ul;
            }
            sum_chi2[run_id] += chi2;
            sum_chi2sq[run_id] += chi2 * chi2;
            ntracks[run_id] += 1ul;
        }
        // if (h_clippingness.find(run_id) == h_clippingness.end()) {
        //     h_clippingness[run_id] = new TH1D(Form("h_clippingness_%d", run_id), Form("Clippingness for run %d; Clippingness; Entries;", run_id), 100, 0.0, 20000.0);
        // }
        if (prvts.find(run_id) == prvts.end()) {
            prvts[run_id] = TTimeStamp{sec, nsec};
            h_time_to_previous_muon[run_id] = new TH1D(Form("h_time_to_previous_muon_%d", run_id), Form("Time to previous muon for run %d; #Delta t (s); Entries;", run_id), 100, 0.0, 5.0);
            continue;
        }
        TTimeStamp ts{sec, nsec};
        h_time_to_previous_muon[run_id]->Fill(ts - prvts[run_id]);
        prvts[run_id] = ts;
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

    for (const auto& [run_id, h] : h_time_to_previous_muon) {
        fit_res[run_id] = fit_rate(h);
    }
    for (const auto& [run_id, h] : h_time_to_previous_muon_single) {
        fit_res_single[run_id] = fit_rate(h);
    }
    for (const auto& [run_id, h] : h_time_to_previous_muon_bundle) {
        fit_res_bundle[run_id] = fit_rate(h);
    }

    TH1D* h_rate_per_run = new TH1D("h_rate_per_run", "Rate per run;RUN ID;Rate (cps);", max_run_id - min_run_id + 1, min_run_id, max_run_id + 1);
    TH1D* h_rate_per_run_single = new TH1D("h_rate_per_run_single", "Rate per run for single;RUN ID;Rate (cps);", max_run_id - min_run_id + 1, min_run_id, max_run_id + 1);
    TH1D* h_rate_per_run_bundle = new TH1D("h_rate_per_run_bundle", "Rate per run for bundle;RUN ID;Rate (cps);", max_run_id - min_run_id + 1, min_run_id, max_run_id + 1);
    for (const auto& [run_id, h] : h_time_to_previous_muon) {
        h_rate_per_run->SetBinContent(run_id - min_run_id + 1, fit_res[run_id]->GetParameter(1));
        h_rate_per_run->SetBinError(run_id - min_run_id + 1, fit_res[run_id]->GetParError(1));
    }
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
    h_rate_per_run_single->SetStats(0);
    h_rate_per_run_single->GetYaxis()->SetMaxDigits(3);
    h_rate_per_run_single->GetXaxis()->CenterTitle(true);
    h_rate_per_run_single->GetYaxis()->CenterTitle(true);
    h_rate_per_run_single->GetYaxis()->SetTitleOffset(1.25);
    h_rate_per_run_single->SetMarkerStyle(kFullCircle);
    h_rate_per_run_single->SetMarkerColor(kViolet);
    h_rate_per_run_single->SetMarkerSize(1.0);
    h_rate_per_run_single->SetLineWidth(2);
    h_rate_per_run_single->SetLineColor(kViolet);
    h_rate_per_run_single->SetMinimum(0.0);
    h_rate_per_run_single->Draw("SAME");
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

    TCanvas* c_chi2_per_run = new TCanvas("c_chi2_per_run", "Chi2 per run", 1000, 1000);
    c_chi2_per_run->cd();
    TH1D* h_chi2_per_run = new TH1D("h_chi2_per_run", "Chi2 per run;RUN ID;Chi2;", max_run_id - min_run_id + 1, min_run_id, max_run_id + 1);
    for (const auto& [run_id, h] : h_time_to_previous_muon) {
        if (ntracks[run_id] == 0) continue;
        h_chi2_per_run->SetBinContent(run_id - min_run_id + 1, sum_chi2[run_id] / ntracks[run_id]);
        // h_chi2_per_run->SetBinError(run_id - min_run_id + 1, sqrt(sum_chi2sq[run_id] / ntracks[run_id] - std::pow(sum_chi2[run_id] / ntracks[run_id], 2.0)));
    }
    h_chi2_per_run->SetStats(0);
    h_chi2_per_run->GetYaxis()->SetMaxDigits(3);
    h_chi2_per_run->GetXaxis()->CenterTitle(true);
    h_chi2_per_run->GetYaxis()->CenterTitle(true);
    h_chi2_per_run->GetYaxis()->SetTitleOffset(1.25);
    h_chi2_per_run->SetMarkerStyle(kFullCircle);
    h_chi2_per_run->SetMarkerColor(kBlue);
    h_chi2_per_run->SetMarkerSize(1.0);
    h_chi2_per_run->SetLineWidth(2);
    h_chi2_per_run->SetLineColor(kBlue);
    h_chi2_per_run->SetMinimum(0.0);
    h_chi2_per_run->Draw();
    c_chi2_per_run->Update();
    c_chi2_per_run->SetTickx();
    c_chi2_per_run->SetTicky();
    c_chi2_per_run->SetGrid();
    c_chi2_per_run->Update();

    // if (int res = fit_and_plot_rate(h_time_to_previous_muon_cdwpttchi2)) return res;

    return 0;
}