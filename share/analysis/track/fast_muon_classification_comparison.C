#include <cstdint>
#include <iostream>

#include <TFile.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TVector3.h>

struct classification {

    int run_id;
    TTimeStamp ts;
    double totq_cd;
    double totq_wp;
    std::size_t ntracks;

};

bool operator<(const classification& a, const classification& b) {
    return a.ts < b.ts;
}

void plot_metrics(const std::map<std::string, TH1D*>& hists, const std::map<std::string, double>& quantiles) {
    if (hists.empty()) return;
    std::map<std::string, Color_t> colors = {
        {"CdWpTtChi2", kBlack}, 
        {"CdClassify", kGreen + 2}, 
        {"WpBasic", kViolet}, 
        {"Amber_v5.5", kBlue},
        {"Edwin", kRed},
    };
    TCanvas* c = new TCanvas(Form("c_%s", hists.at("CdWpTtChi2")->GetName()), "Metric", 1000, 1000);
    c->cd();

    double max = 0.0;
    for (const auto& [method, h] : hists) {
        if (h->GetMaximum() > max) {
            max = h->GetMaximum();
        }
        std::cout << "68.2% " << h->GetName() << " (" << h->GetEntries() << " entries): " << quantiles.at(method) << '\n';
    }

    TLegend* leg = new TLegend(0.45, 0.65, 0.85, 0.85);

    bool first = true;
    for (const auto& [method, h] : hists) {
        h->SetStats(0);
        h->SetLineColor(colors[method]);
        h->SetLineWidth(3);
        h->GetXaxis()->SetMaxDigits(3);
        h->GetYaxis()->SetMaxDigits(3);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->GetYaxis()->SetTitleOffset(1.25);
        if (first) {
            first = false;
            h->SetMaximum(max * 1.1);
            h->Draw();
        }
        else {
            h->Draw("SAME");
        }

        TLine* line = new TLine(quantiles.at(method), 0.0, quantiles.at(method), max * 1.1);
        line->SetLineStyle(2);
        line->SetLineWidth(3);
        line->SetLineColor(colors[method]);
        line->Draw();

        leg->AddEntry(h, Form("%s: 68%% quantile = %.2f", method.c_str(), quantiles.at(method)), "l");
    }
    leg->SetTextSize(0.02);
    leg->Draw();

    c->SetTickx();
    c->SetTicky();
    c->SetGrid();
    c->Update();
}

std::set<classification> open_amber_v5_5_user_chain(const char* path) {
    TChain* chain = new TChain("MuonReco");
    chain->Add(path);
    std::set<track> classifications;

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

    long nentries = chain->GetEntries();
    std::cout << "Info: Found " << nentries << " entries in Amber_v5.5 files\n";
    for (long k = 0l; k < nentries; ++k) {
        chain->GetEntry(k);
        classifications.insert(classification{
            .run_id = runID,
            .ts = TTimeStamp(fSec, fNanoSec),
            .totq_cd = 0.0,
            .totq_wp = charge,
            .ntracks = (muonType == 0 ? 1ul : 2ul)
        });
    }
    return classifications;
}

std::map<std::string, std::set<classification>> open_joint_reco_user_chain(const char* path) {
    TChain* chain = new TChain("muons");
    chain->Add(path);
    std::map<std::string, std::set<classification>> classifications;

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

    long nentries = chain->GetEntries();
    std::cout << "Info: Found " << nentries << " entries in joint reco files\n";
    for (long k = 0l; k < nentries; ++k) {
        chain->GetEntry(k);
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
        for (std::size_t i = 0ul; i < method->size(); ++i) {
            if ((*method)[i] != "CdClassify" || (*method)[i] != "WpBasic") continue;
            classifications[(*method)[i]].insert(classification{
                .run_id = run_id,
                .ts = TTimeStamp(sec, nsec),
                .totq_cd = totq_cd,
                .totq_wp = totq_wp,
                .ntracks = (*method)[i] == "CdClassify" ? ntracks_cdclassify : ntracks_wpclassify
            });
        }
    }
    return classifications;
}


struct MuonClassification {
    int run_id;
    std::size_t ntracks_wpclassify;
    std::size_t ntracks_cdclassify;
    std::size_t ntracks_amber;
};

std::vector<MuonClassification> compute_global_correlations_classification(std::map<std::string, std::set<classification>>& classifications) {
    if (classifications.find("WpBasic") == classifications.end()) {
        std::cerr << "Error: WpBasic classifications not found in map\n";
        return {};
    }
    const std::set<classification>& wp_classifications = classifications["WpBasic"];

    std::vector<MuonClassification> results;

    for (const classification& wp_clas : wp_classifications) {
        std::map<std::string, classification> coincident_map;
        std::size_t ntracks_cdclassify = 0ul;
        std::size_t ntracks_amber = 0ul;

        for (const auto& [method, clas_set] : classifications) {
            if (method == "WpBasic") continue;
            bool found_in_method = false;
            TTimeStamp lower_bound_ts(wp_clas.ts.GetSec(), wp_clas.ts.GetNanoSec() - 1000);
            TTimeStamp upper_bound_ts(wp_clas.ts.GetSec(), wp_clas.ts.GetNanoSec() + 1000);
            std::set<classification>::const_iterator it = clas_set.lower_bound({0, lower_bound_ts, 0, 0, {}, {}});
            
            while (it != clas_set.end() && lower_bound_ts <= it->ts && it->ts <= upper_bound_ts) {
                coincident_map[method] = *it;
                found_in_method = true;
                if (method == "CdClassify") {
                    ntracks_cdclassify = it->ntracks;
                }
                if (method == "Amber_v5.5") {
                    ++ntracks_amber = it->ntracks;
                }
                ++it;
            }

            if (!found_in_method) {
                all_found = false;
                break; 
            }
        }

        if (all_found) {
            classifications.push_back(MuonClassification{
                .run_id = wp_clas.run_id,
                .ntracks_wpclassify = wp_clas.ntracks,
                .ntracks_cdclassify = ntracks_cdclassify,
                .ntracks_amber = ntracks_amber
            });
        }
    }
    return classifications;
}

int fast_muon_reconstruction_comparison(const char* path_joint, const char* path_cdwpttchi2, const char* path_amber, const char* path_edwin) {
    std::map<std::string, std::set<classification>> classifications = open_joint_reco_user_chain(path_joint);
    classifications["Amber_v5.5"] = open_amber_v5_5_user_chain(path_amber);

    std::vector<MuonClassification> results = compute_global_correlations_classification(classifications);

    std::cout << "Results size: " << results.size() << '\n';

    return 0;
}