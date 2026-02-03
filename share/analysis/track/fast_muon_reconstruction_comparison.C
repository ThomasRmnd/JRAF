#include <iostream>

#include <TFile.h>
#include <TTimeStamp.h>
#include <TTree.h>

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

    TH1D* h_angle = new TH1D("h_angle", "Angle between tracks direction;#alpha (deg);Entries;", 25, 0.0, 5.0);
    TH1D* h_distance = new TH1D("h_distance", "Distance between tracks middle point;d_{mid} (m);Entries;", 25, 0.0, 2.0);

    long nentries = tree->GetEntries();
    for (long k = 0l; k < nentries; ++k) {
        tree->GetEntry(k);
        bool found_cdwptt = false, found_tt = false;
        std::size_t k_cdwptt = 0ul, k_tt = 0ul;
        std::size_t n_cdclassify = 0ul;
        for (std::size_t i = 0ul; i < method->size(); ++i) {
            if ((*method)[i] == "CdWpTtChi2") {
                k_cdwptt = i;
                found_cdwptt = true;
            }
            else if ((*method)[i] == "Tt") {
                k_tt = i;
                found_tt = true;
            }
            else if ((*method)[i] == "CdClassify") {
                ++n_cdclassify;
            }
        }
        if (!found_cdwptt || !found_tt) continue;
        if (n_cdclassify != 1) continue;
        TVector3 ipos_cdwptt((*iposx)[k_cdwptt], (*iposy)[k_cdwptt], (*iposz)[k_cdwptt]);
        TVector3 ipos_tt((*iposx)[k_tt], (*iposy)[k_tt], (*iposz)[k_tt]);
        TVector3 fpos_cdwptt((*fposx)[k_cdwptt], (*fposy)[k_cdwptt], (*fposz)[k_cdwptt]);
        TVector3 fpos_tt((*fposx)[k_tt], (*fposy)[k_tt], (*fposz)[k_tt]);
        TVector3 dir_cdwptt = (fpos_cdwptt - ipos_cdwptt).Unit();
        TVector3 dir_tt = (fpos_tt - ipos_tt).Unit();
        TVector3 mpos_cdwptt = (ipos_cdwptt + fpos_cdwptt) * 0.5;
        TVector3 mpos_tt = (ipos_tt + fpos_tt) * 0.5;
        double angle = dir_cdwptt.Angle(dir_tt);
        double distance = (mpos_cdwptt - mpos_tt).Mag();
        h_angle->Fill(angle * 180.0 / TMath::Pi());
        h_distance->Fill(distance / 1000.0);
    }

    TCanvas* c_angle = new TCanvas("c_angle", "Angle between tracks direction", 1000, 1000);
    c_angle->cd();

    h_angle->SetStats(0);
    h_angle->SetLineColor(kBlue);
    h_angle->SetLineWidth(3);
    h_angle->GetXaxis()->SetMaxDigits(3);
    h_angle->GetYaxis()->SetMaxDigits(3);
    h_angle->GetXaxis()->CenterTitle(true);
    h_angle->GetYaxis()->CenterTitle(true);
    h_angle->GetYaxis()->SetTitleOffset(1.25);
    h_angle->Draw();

    c_angle->SetTickx();
    c_angle->SetTicky();
    c_angle->SetGrid();
    c_angle->Update();

    TCanvas* c_distance = new TCanvas("c_distance", "Distance between tracks middle point", 1000, 1000);
    c_distance->cd();

    h_distance->SetStats(0);
    h_distance->SetLineColor(kBlue);
    h_distance->SetLineWidth(3);
    h_distance->GetXaxis()->SetMaxDigits(3);
    h_distance->GetYaxis()->SetMaxDigits(3);
    h_distance->GetXaxis()->CenterTitle(true);
    h_distance->GetYaxis()->CenterTitle(true);
    h_distance->GetYaxis()->SetTitleOffset(1.25);
    h_distance->Draw();

    c_distance->SetTickx();
    c_distance->SetTicky();
    c_distance->SetGrid();
    c_distance->Update();

    return 0;
}