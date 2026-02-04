#include <iostream>

#include <TFile.h>
#include <TTimeStamp.h>
#include <TTree.h>

double get_quantile(std::vector<double>::const_iterator first, std::vector<double>::const_iterator last, double quantile) {
    if (first == last) return 0.0;
    std::vector<double> v(first, last);
    std::vector<double>::iterator position = v.begin() + quantile * v.size();
    std::nth_element(v.begin(), position, v.end());
    return *position;
}

void plot_metrics(std::vector<TH1D*> hists) {
    if (hists.empty()) return;
    std::vector<Color_t> colors = {kBlack, kGreen + 2, kViolet};
    TCanvas* c = new TCanvas(Form("c_%s", hists[0]->GetName()), "Metric", 1000, 1000);
    c->cd();

    for (std::size_t i = 0ul; i < hists.size(); ++i) {
        TH1D* h = hists[i];
        
        h->SetStats(0);
        h->SetLineColor(colors[i]);
        h->SetLineWidth(3);
        h->GetXaxis()->SetMaxDigits(3);
        h->GetYaxis()->SetMaxDigits(3);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->GetYaxis()->SetTitleOffset(1.25);
        if (i == 0) {
            h->Draw();
        }
        else {
            h->Draw("SAME");
        }
    };

    c->SetTickx();
    c->SetTicky();
    c->SetGrid();
    c->Update();
}

int fast_muon_reconstruction_comparison(const char* filepath) {
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
    std::vector<std::string> *method = nullptr;
    std::vector<unsigned char> *det = nullptr;
    std::vector<double> *quality = nullptr;
    std::vector<double> *iposx = nullptr, *iposy = nullptr, *iposz = nullptr; 
    std::vector<double> *fposx = nullptr, *fposy = nullptr, *fposz = nullptr;
    tree->SetBranchAddress("run_id", &run_id);
    tree->SetBranchAddress("sec", &sec);
    tree->SetBranchAddress("nsec", &nsec);
    // tree->SetBranchAddress("totq_cd", &totq_cd);
    // tree->SetBranchAddress("totq_wp", &totq_wp);
    tree->SetBranchAddress("method", &method);
    tree->SetBranchAddress("det", &det);
    tree->SetBranchAddress("quality", &quality);
    tree->SetBranchAddress("iposx", &iposx);
    tree->SetBranchAddress("iposy", &iposy);
    tree->SetBranchAddress("iposz", &iposz);
    tree->SetBranchAddress("fposx", &fposx);
    tree->SetBranchAddress("fposy", &fposy);
    tree->SetBranchAddress("fposz", &fposz);

    std::map<std::string, TH1D*> method_angle_map;
    std::map<std::string, TH1D*> method_distance_map;

    method_angle_map["CdWpTtChi2"] = new TH1D("h_angle_cdwpttchi2", "Angle between tracks direction (CdWpTtChi2);#alpha (deg);Entries;", 25, 0.0, 5.0);
    method_distance_map["CdWpTtChi2"] = new TH1D("h_distance_cdwpttchi2", "Distance between tracks middle point (CdWpTtChi2);d_{mid} (m);Entries;", 25, 0.0, 2.0);
    method_angle_map["CdClassify"] = new TH1D("h_angle_cdclassify", "Angle between tracks direction (CdClassify);#alpha (deg);Entries;", 25, 0.0, 5.0);
    method_distance_map["CdClassify"] = new TH1D("h_distance_cdclassify", "Distance between tracks middle point (CdClassify);d_{mid} (m);Entries;", 25, 0.0, 2.0);
    method_angle_map["WpBasic"] = new TH1D("h_angle_wpclassify", "Angle between tracks direction (WpClassify);#alpha (deg);Entries;", 25, 0.0, 5.0);
    method_distance_map["WpBasic"] = new TH1D("h_distance_wpclassify", "Distance between tracks middle point (WpClassify);d_{mid} (m);Entries;", 25, 0.0, 2.0);

    std::vector<double> angles;
    std::vector<double> distances;

    long nentries = tree->GetEntries();
    for (long k = 0l; k < nentries; ++k) {
        tree->GetEntry(k);
        std::map<std::string, std::vector<std::size_t>> track_method_map;
        for (std::size_t i = 0ul; i < method->size(); ++i) {
            track_method_map[(*method)[i]].push_back(i);
        }
        if (track_method_map.find("CdWpTtChi2") == track_method_map.end()) continue;
        if (track_method_map.find("CdClassify") == track_method_map.end()) continue;
        if (track_method_map.find("Tt") == track_method_map.end()) continue;

        if (track_method_map["Tt"].size() != 1) continue;
        if (track_method_map["CdClassify"].size() != 1) continue;

        std::size_t k_tt = track_method_map["Tt"][0];
        TVector3 ipos_tt((*iposx)[k_tt], (*iposy)[k_tt], (*iposz)[k_tt]);
        TVector3 fpos_tt((*fposx)[k_tt], (*fposy)[k_tt], (*fposz)[k_tt]);
        TVector3 dir_tt = (fpos_tt - ipos_tt).Unit();
        TVector3 mpos_tt = (ipos_tt + fpos_tt) * 0.5;

        if (mpos_tt.Mag() > 17700.0) continue;

        for (const auto& [key, val] : track_method_map) {
            if (key == "Tt") continue;
            for (std::vector<std::size_t>::const_iterator it = val.begin(); it != val.end(); ++it) {
                TVector3 ipos((*iposx)[*it], (*iposy)[*it], (*iposz)[*it]);
                TVector3 fpos((*fposx)[*it], (*fposy)[*it], (*fposz)[*it]);
                TVector3 dir = (fpos - ipos).Unit();
                TVector3 mpos = (ipos + fpos) * 0.5;

                double angle = dir.Angle(dir_tt) * 180.0 / TMath::Pi();
                double distance = (mpos - mpos_tt).Mag() / 1000.0;
                method_angle_map[key]->Fill(angle);
                method_distance_map[key]->Fill(distance);
                if (key == "CdWpTtChi2") {
                    angles.push_back(angle);
                    distances.push_back(distance);
                }
            }
        }
    }

    std::cout << "68.2% angle: " << get_quantile(angles.begin(), angles.end(), 0.682) << '\n';
    std::cout << "95.4% angle: " << get_quantile(angles.begin(), angles.end(), 0.954) << '\n';
    std::cout << "68.2% distance: " << get_quantile(distances.begin(), distances.end(), 0.682) << '\n';
    std::cout << "95.4% distance: " << get_quantile(distances.begin(), distances.end(), 0.954) << '\n';

    plot_metrics({method_angle_map["CdWpTtChi2"], method_angle_map["CdClassify"], method_angle_map["WpBasic"]});
    plot_metrics({method_distance_map["CdWpTtChi2"], method_distance_map["CdClassify"], method_distance_map["WpBasic"]});

    return 0;
}