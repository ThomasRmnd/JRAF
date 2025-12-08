#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <unordered_set>

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TVector3.h>

enum class StatOpt {
    None        =         0,
    Name        =         1,  // n = 1
    Entries     =        10,  // e = 1
    Mean        =       100,  // m = 1
    MeanError   =       200,  // m = 2
    RMS         =      1000,  // r = 1
    RMSError    =      2000,  // r = 2
    Underflow   =     10000,  // u = 1
    Overflow    =    100000,  // o = 1
    Integral    =   1000000,  // i = 1
    Skew        =  10000000,  // s = 1
    SkewError   =  20000000,  // s = 2
    Kurtosis    = 100000000,  // k = 1
    KurtError   = 200000000,  // k = 2
};

inline StatOpt operator|(StatOpt a, StatOpt b) {
    return static_cast<StatOpt>(static_cast<int>(a) + static_cast<int>(b));
}

inline int ToROOTOpt(StatOpt mode) {
    return static_cast<int>(mode);
}

enum class FitOpt {
    None        =    0,
    Params      =    1,  // v = 1 (requires FitErrors)
    AllParams   =    2,  // v = 2
    Errors      =   10,  // e = 1
    Chi2NDF     =  100,  // c = 1
    Proba       = 1000,  // p = 1
};

inline FitOpt operator|(FitOpt a, FitOpt b) {
    return static_cast<FitOpt>(static_cast<int>(a) + static_cast<int>(b));
}

inline int ToROOTOpt(FitOpt mode) {
    return static_cast<int>(mode);
}

inline TTimeStamp operator-(const TTimeStamp& lhs, const TTimeStamp& rhs) {
    return TTimeStamp(lhs.GetSec() - rhs.GetSec(), lhs.GetNanoSec() - rhs.GetNanoSec());
}

inline TTimeStamp operator+(const TTimeStamp& lhs, const TTimeStamp& rhs) {
    return TTimeStamp(lhs.GetSec() + rhs.GetSec(), lhs.GetNanoSec() + rhs.GetNanoSec());
}

void no_good_tt_maches(const char* filename) {
    TFile* file = TFile::Open(filename, "READ");
    if (!file) {
        std::cerr << "Cannot open file " << filename << '\n';
        return;
    }
    TTree* tree = file->Get<TTree>("matches");
    if (!tree) {
        std::cerr << "Cannot retrieve tree matches in file " << filename << '\n';
        return;
    }

    Int_t run_id;
    Long_t sec;
    Int_t nsec;
    std::vector<std::string>* method = nullptr;
    std::vector<unsigned char>* det = nullptr;
    std::vector<double>* quality = nullptr;
    std::vector<double> *iposx = nullptr, *iposy = nullptr, *iposz = nullptr, *fposx = nullptr, *fposy = nullptr, *fposz = nullptr;

    tree->SetBranchAddress("run_id", &run_id);
    tree->SetBranchAddress("sec", &sec);
    tree->SetBranchAddress("nsec", &nsec);
    tree->SetBranchAddress("method", &method);
    tree->SetBranchAddress("det", &det);
    tree->SetBranchAddress("quality", &quality);
    tree->SetBranchAddress("iposx", &iposx);
    tree->SetBranchAddress("iposy", &iposy);
    tree->SetBranchAddress("iposz", &iposz);
    tree->SetBranchAddress("fposx", &fposx);
    tree->SetBranchAddress("fposy", &fposy);
    tree->SetBranchAddress("fposz", &fposz);

    Int_t evtID, NTotPoints, NTracks;
    Int_t NPoints[20];
    TTimeStamp* start_TS = nullptr;
    Float_t PointX[20], PointY[20], PointZ[20];
    Double_t Coeff0[20], Coeff1[20], Coeff2[20], Coeff3[20], Coeff4[20], Coeff5[20];
    Double_t Chi2[20];

    tree->SetBranchAddress("evtID", &evtID);
    tree->SetBranchAddress("NTotPoints", &NTotPoints);
    tree->SetBranchAddress("PointX", &PointX);
    tree->SetBranchAddress("PointY", &PointY);
    tree->SetBranchAddress("PointZ", &PointZ);
    tree->SetBranchAddress("NTracks", &NTracks);
    tree->SetBranchAddress("NPoints", NPoints);
    tree->SetBranchAddress("start_TS", &start_TS);
    tree->SetBranchAddress("Coeff0", &Coeff0);
    tree->SetBranchAddress("Coeff1", &Coeff1);
    tree->SetBranchAddress("Coeff2", &Coeff2);
    tree->SetBranchAddress("Coeff3", &Coeff3);
    tree->SetBranchAddress("Coeff4", &Coeff4);
    tree->SetBranchAddress("Coeff5", &Coeff5);
    tree->SetBranchAddress("Chi2", &Chi2);

    TVector3 ipos, fpos;
    TVector3 pos_cdwp, dir_cdwp, pos_tt, dir_tt;

    std::unordered_set<int> run_with_tt;
    for (int k = 0; k < tree->GetEntries(); ++k) {
        if (k % 1000 == 0) std::cout << "\rEntries: " << k << " / " << tree->GetEntries();
        tree->GetEntry(k);
        if (NTotPoints == 0) continue;
        run_with_tt.insert(run_id);
    }

    std::cout << "[INFO] Number of entries: " << tree->GetEntries() << '\n';
    std::size_t nb_cdwp_reco = 0ul;
    bool set_cdwp_previous_muon = false;
    TTimeStamp cdwp_previous_muon{0, 0};
    TH1D* h_cdwp_muon_rate = new TH1D("h_cdwp_muon_rate", "h_cdwp_muon_rate", 50, 0.0, 5.0);
    bool set_cdwp_tt_3pts_previous_muon = false;
    bool set_cdwp_tt_2pts_previous_muon = false;
    TTimeStamp cdwp_tt_3pts_previous_muon{0, 0};
    TTimeStamp cdwp_tt_2pts_previous_muon{0, 0};
    TH1D* h_cdwp_tt_3pts_muon_rate = new TH1D("h_cdwp_tt_3pts_muon_rate", "h_cdwp_tt_3pts_muon_rate", 50, 0.0, 20.0);
    TH1D* h_cdwp_tt_2pts_muon_rate = new TH1D("h_cdwp_tt_2pts_muon_rate", "h_cdwp_tt_2pts_muon_rate", 50, 0.0, 20.0);
    bool set_cdwp_no_tt_rate_previous_muon = false;
    TTimeStamp cdwp_no_tt_rate_previous_muon{0, 0};
    TH1D* h_cdwp_no_tt_rate = new TH1D("h_cdwp_no_tt_rate", "h_cdwp_no_tt_rate", 50, 0.0, 20.0);
    for (int k = 0; k < tree->GetEntries(); ++k) {
        if (k % 1000 == 0) std::cout << "\rEntries: " << k << " / " << tree->GetEntries();
        tree->GetEntry(k);

        if (run_with_tt.find(run_id) == run_with_tt.end()) continue;

        TTimeStamp ts{sec, nsec};

        std::unordered_map<std::string, std::size_t> method_count = {
            {"CdWpTtChi2", 0ul},
            {"WpClassify", 0ul}
        };
        std::size_t j = 0ul;
        for (std::size_t i = 0ul; i < method->size(); ++i) {
            ++(method_count[(*method)[i]]);
            if ((*method)[i] == "CdWpTtChi2") j = i;
        }
        if (/* method_count["WpClassify"] != 1ul || */ method_count["CdWpTtChi2"] != 1ul) continue;

        ipos.SetXYZ((*iposx)[j], (*iposy)[j], (*iposz)[j]);
        fpos.SetXYZ((*fposx)[j], (*fposy)[j], (*fposz)[j]);

        pos_cdwp = ipos;
        dir_cdwp = (fpos - ipos).Unit();

        pos_tt.SetXYZ(Coeff0[0], Coeff1[0], Coeff2[0] + 26452.0);
        dir_tt.SetXYZ(Coeff3[0], Coeff4[0], Coeff5[0]);
        dir_tt = dir_tt.Unit();

        if (!set_cdwp_previous_muon) {
            set_cdwp_previous_muon = true;
            cdwp_previous_muon = ts;
        }
        else {
            h_cdwp_muon_rate->Fill(ts - cdwp_previous_muon);
            cdwp_previous_muon = ts;
        }
        ++nb_cdwp_reco;
        if (NTotPoints < 2 /* == 0 */) {
            double layer_tmp[3] = {24500.0, 26000.0, 27500.0};
            double cdwp_reco_x_at_layer[3];
            double cdwp_reco_y_at_layer[3];
            for (int i = 0; i < 3; ++i) {
                double t_cdwp = (layer_tmp[i] - pos_cdwp.Z()) / dir_cdwp.Z();
                TVector3 p_cdwp_at_z = pos_cdwp + t_cdwp * dir_cdwp;
                cdwp_reco_x_at_layer[i] = p_cdwp_at_z.X();
                cdwp_reco_y_at_layer[i] = p_cdwp_at_z.Y();
            }
            if (
                -23500.0 <= cdwp_reco_x_at_layer[0] && cdwp_reco_x_at_layer[0] <= 23500.0 &&
                -10000.0 <= cdwp_reco_y_at_layer[0] && cdwp_reco_y_at_layer[0] <= 10000.0 &&
                -23500.0 <= cdwp_reco_x_at_layer[1] && cdwp_reco_x_at_layer[1] <= 23500.0 &&
                -10000.0 <= cdwp_reco_y_at_layer[1] && cdwp_reco_y_at_layer[1] <= 10000.0 &&
                -23500.0 <= cdwp_reco_x_at_layer[2] && cdwp_reco_x_at_layer[2] <= 23500.0 &&
                -10000.0 <= cdwp_reco_y_at_layer[2] && cdwp_reco_y_at_layer[2] <= 10000.0
            ) {
                if (!set_cdwp_no_tt_rate_previous_muon) {
                    set_cdwp_no_tt_rate_previous_muon = true;
                    cdwp_no_tt_rate_previous_muon = ts;
                }
                else {
                    h_cdwp_no_tt_rate->Fill(ts - cdwp_no_tt_rate_previous_muon);
                    cdwp_no_tt_rate_previous_muon = ts;
                }
            }
        }
        if (NTracks != 1) continue;
        if (NPoints[0] >= 3) {
            if (!set_cdwp_tt_3pts_previous_muon) {
                set_cdwp_tt_3pts_previous_muon = true;
                cdwp_tt_3pts_previous_muon = ts;
            }
            else {
                h_cdwp_tt_3pts_muon_rate->Fill(ts - cdwp_tt_3pts_previous_muon);
                cdwp_tt_3pts_previous_muon = ts;
            }
        }
        else {
            if (!set_cdwp_tt_2pts_previous_muon) {
                set_cdwp_tt_2pts_previous_muon = true;
                cdwp_tt_2pts_previous_muon = ts;
            }
            else {
                h_cdwp_tt_2pts_muon_rate->Fill(ts - cdwp_tt_2pts_previous_muon);
                cdwp_tt_2pts_previous_muon = ts;
            }
        }
    }
    std::cout << "[INFO] Number of CD/WP reco entries: " << nb_cdwp_reco << '\n';

    // ============================================================================================
    // dt previous muon
    // ============================================================================================

    TCanvas* c_cdwp_muon_rate = new TCanvas("c_cdwp_muon_rate", "c_cdwp_muon_rate", 1000, 1000);
    c_cdwp_muon_rate->cd();

    TF1* f_cdwp_exp = new TF1("f_cdwp_exp", "[0] * exp(- [1] * x)", 0.0, 5.0);
    f_cdwp_exp->SetParameter(0, h_cdwp_muon_rate->GetMaximum());
    f_cdwp_exp->SetParameter(1, 5.0);

    f_cdwp_exp->SetLineColor(kRed);
    f_cdwp_exp->SetLineWidth(3);

    h_cdwp_muon_rate->Fit(f_cdwp_exp, "R");
    h_cdwp_muon_rate->SetLineWidth(3);
    h_cdwp_muon_rate->GetXaxis()->SetTitle("#Delta t_{muon} (s)");
    h_cdwp_muon_rate->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_muon_rate->GetYaxis()->SetTitle("Entries");
    h_cdwp_muon_rate->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_muon_rate->GetYaxis()->SetTitleOffset(1.5);
    h_cdwp_muon_rate->Draw();
    f_cdwp_exp->Draw("SAME");
    
    c_cdwp_muon_rate->SetTickx();
    c_cdwp_muon_rate->SetTicky();
    c_cdwp_muon_rate->SetLogy();

    c_cdwp_muon_rate->Update();

    TPaveStats* st_cdwp_muon_rate = (TPaveStats*)h_cdwp_muon_rate->FindObject("stats");
    st_cdwp_muon_rate->SetOptStat(ToROOTOpt(StatOpt::Entries | StatOpt::Mean | StatOpt::RMS));
    st_cdwp_muon_rate->SetOptFit(ToROOTOpt(FitOpt::Proba | FitOpt::Chi2NDF | FitOpt::AllParams | FitOpt::Errors));
    st_cdwp_muon_rate->SetX1NDC(0.5);
    st_cdwp_muon_rate->SetX2NDC(0.85);
    st_cdwp_muon_rate->SetY1NDC(0.5);
    st_cdwp_muon_rate->SetY2NDC(0.85);

    c_cdwp_muon_rate->Modified();
    c_cdwp_muon_rate->Update();

    // ============================================================================================
    // dt previous muon with TT 3 points
    // ============================================================================================

    TCanvas* c_cdwp_tt_3pts_muon_rate = new TCanvas("c_cdwp_tt_3pts_muon_rate", "c_cdwp_tt_3pts_muon_rate", 1000, 1000);
    c_cdwp_tt_3pts_muon_rate->cd();

    TF1* f_cdwp_tt_3pts_exp = new TF1("f_cdwp_tt_3pts_exp", "[0] * exp(- [1] * x)", 0.0, 20.0);
    f_cdwp_tt_3pts_exp->SetParameter(0, h_cdwp_tt_3pts_muon_rate->GetMaximum());
    f_cdwp_tt_3pts_exp->SetParameter(1, 5.0);

    f_cdwp_tt_3pts_exp->SetLineColor(kRed);
    f_cdwp_tt_3pts_exp->SetLineWidth(3);

    h_cdwp_tt_3pts_muon_rate->Fit(f_cdwp_tt_3pts_exp, "R");
    h_cdwp_tt_3pts_muon_rate->SetLineWidth(3);
    h_cdwp_tt_3pts_muon_rate->GetXaxis()->SetTitle("#Delta t_{muon} (s)");
    h_cdwp_tt_3pts_muon_rate->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_3pts_muon_rate->GetYaxis()->SetTitle("Entries");
    h_cdwp_tt_3pts_muon_rate->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_3pts_muon_rate->GetYaxis()->SetTitleOffset(1.5);
    h_cdwp_tt_3pts_muon_rate->Draw();
    f_cdwp_tt_3pts_exp->Draw("SAME");
    
    c_cdwp_tt_3pts_muon_rate->SetTickx();
    c_cdwp_tt_3pts_muon_rate->SetTicky();
    c_cdwp_tt_3pts_muon_rate->SetLogy();

    c_cdwp_tt_3pts_muon_rate->Update();

    TPaveStats* st_cdwp_tt_3pts_muon_rate = (TPaveStats*)h_cdwp_tt_3pts_muon_rate->FindObject("stats");
    st_cdwp_tt_3pts_muon_rate->SetOptStat(ToROOTOpt(StatOpt::Entries | StatOpt::Mean | StatOpt::RMS));
    st_cdwp_tt_3pts_muon_rate->SetOptFit(ToROOTOpt(FitOpt::Proba | FitOpt::Chi2NDF | FitOpt::AllParams | FitOpt::Errors));
    st_cdwp_tt_3pts_muon_rate->SetX1NDC(0.5);
    st_cdwp_tt_3pts_muon_rate->SetX2NDC(0.85);
    st_cdwp_tt_3pts_muon_rate->SetY1NDC(0.5);
    st_cdwp_tt_3pts_muon_rate->SetY2NDC(0.85);

    c_cdwp_tt_3pts_muon_rate->Modified();
    c_cdwp_tt_3pts_muon_rate->Update();

    // ============================================================================================
    // dt previous muon with TT 2 points
    // ============================================================================================

    TCanvas* c_cdwp_tt_2pts_muon_rate = new TCanvas("c_cdwp_tt_2pts_muon_rate", "c_cdwp_tt_2pts_muon_rate", 1000, 1000);
    c_cdwp_tt_2pts_muon_rate->cd();

    TF1* f_cdwp_tt_2pts_exp = new TF1("f_cdwp_tt_2pts_exp", "[0] * exp(- [1] * x)", 0.0, 20.0);
    f_cdwp_tt_2pts_exp->SetParameter(0, h_cdwp_tt_2pts_muon_rate->GetMaximum());
    f_cdwp_tt_2pts_exp->SetParameter(1, 5.0);

    f_cdwp_tt_2pts_exp->SetLineColor(kRed);
    f_cdwp_tt_2pts_exp->SetLineWidth(3);

    h_cdwp_tt_2pts_muon_rate->Fit(f_cdwp_tt_2pts_exp, "R");
    h_cdwp_tt_2pts_muon_rate->SetLineWidth(3);
    h_cdwp_tt_2pts_muon_rate->GetXaxis()->SetTitle("#Delta t_{muon} (s)");
    h_cdwp_tt_2pts_muon_rate->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_2pts_muon_rate->GetYaxis()->SetTitle("Entries");
    h_cdwp_tt_2pts_muon_rate->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_2pts_muon_rate->GetYaxis()->SetTitleOffset(1.5);
    h_cdwp_tt_2pts_muon_rate->Draw();
    f_cdwp_tt_2pts_exp->Draw("SAME");
    
    c_cdwp_tt_2pts_muon_rate->SetTickx();
    c_cdwp_tt_2pts_muon_rate->SetTicky();
    c_cdwp_tt_2pts_muon_rate->SetLogy();

    c_cdwp_tt_2pts_muon_rate->Update();

    TPaveStats* st_cdwp_tt_2pts_muon_rate = (TPaveStats*)h_cdwp_tt_2pts_muon_rate->FindObject("stats");
    st_cdwp_tt_2pts_muon_rate->SetOptStat(ToROOTOpt(StatOpt::Entries | StatOpt::Mean | StatOpt::RMS));
    st_cdwp_tt_2pts_muon_rate->SetOptFit(ToROOTOpt(FitOpt::Proba | FitOpt::Chi2NDF | FitOpt::AllParams | FitOpt::Errors));
    st_cdwp_tt_2pts_muon_rate->SetX1NDC(0.5);
    st_cdwp_tt_2pts_muon_rate->SetX2NDC(0.85);
    st_cdwp_tt_2pts_muon_rate->SetY1NDC(0.5);
    st_cdwp_tt_2pts_muon_rate->SetY2NDC(0.85);

    c_cdwp_tt_2pts_muon_rate->Modified();
    c_cdwp_tt_2pts_muon_rate->Update();

    // ============================================================================================
    // dt previous muon with no TT
    // ============================================================================================

    TCanvas* c_cdwp_no_tt_rate = new TCanvas("c_cdwp_no_tt_rate", "c_cdwp_no_tt_rate", 1000, 1000);
    c_cdwp_no_tt_rate->cd();

    TF1* f_cdwp_no_tt_exp = new TF1("f_cdwp_no_tt_exp", "[0] * exp(-x / [1]) + [2] * exp(-x / [3])", 0.0, 20.0);
    double maxv = h_cdwp_no_tt_rate->GetMaximum();
    f_cdwp_no_tt_exp->SetParameters(maxv, 1.0, maxv * 0.1, 5.0);

    f_cdwp_no_tt_exp->SetLineColor(kRed);
    f_cdwp_no_tt_exp->SetLineWidth(3);

    h_cdwp_no_tt_rate->Fit(f_cdwp_no_tt_exp, "R");
    h_cdwp_no_tt_rate->SetLineWidth(3);
    h_cdwp_no_tt_rate->GetXaxis()->SetTitle("#Delta t_{muon} (s)");
    h_cdwp_no_tt_rate->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_no_tt_rate->GetYaxis()->SetTitle("Entries");
    h_cdwp_no_tt_rate->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_no_tt_rate->GetYaxis()->SetTitleOffset(1.5);
    h_cdwp_no_tt_rate->Draw();
    f_cdwp_no_tt_exp->Draw("SAME");
    
    c_cdwp_no_tt_rate->SetTickx();
    c_cdwp_no_tt_rate->SetTicky();
    c_cdwp_no_tt_rate->SetLogy();

    c_cdwp_no_tt_rate->Update();

    TPaveStats* st_cdwp_no_tt_rate = (TPaveStats*)h_cdwp_no_tt_rate->FindObject("stats");
    st_cdwp_no_tt_rate->SetOptStat(ToROOTOpt(StatOpt::Entries | StatOpt::Mean | StatOpt::RMS));
    st_cdwp_no_tt_rate->SetOptFit(ToROOTOpt(FitOpt::Proba | FitOpt::Chi2NDF | FitOpt::AllParams | FitOpt::Errors));
    st_cdwp_no_tt_rate->SetX1NDC(0.5);
    st_cdwp_no_tt_rate->SetX2NDC(0.85);
    st_cdwp_no_tt_rate->SetY1NDC(0.5);
    st_cdwp_no_tt_rate->SetY2NDC(0.85);

    c_cdwp_no_tt_rate->Modified();
    c_cdwp_no_tt_rate->Update();

    // ============================================================================================
    // 
    // ============================================================================================

    TH1D* h_ts_diff = new TH1D("h_ts_diff", "h_ts_diff", 200, -1000.0, 1000.0);
    TH1D* h_ts_diff_1trk = new TH1D("h_ts_diff_1trk", "h_ts_diff_1trk", 200, -1000.0, 1000.0);
    TH1D* h_ts_diff_1trk_2pts = new TH1D("h_ts_diff_1trk_2pts", "h_ts_diff_1trk_2pts", 200, -1000.0, 1000.0);
    TH1D* h_ts_diff_1trk_3pts = new TH1D("h_ts_diff_1trk_3pts", "h_ts_diff_1trk_3pts", 200, -1000.0, 1000.0);
    TH1D* h_ts_diff_1trk_3pts_dlay = new TH1D("h_ts_diff_1trk_3pts_dlay", "h_ts_diff_1trk_3pts_dlay", 200, -1000.0, 1000.0);

    int nb_bins_x_tt = 926; // 1850; // 100;
    double xmin_tt = -25.0;
    double xmax_tt = 25.0;
    int nb_bins_y_tt = 444; // 890; // 100
    double ymin_tt = -12.0;
    double ymax_tt = 12.0;
    
    TH1D* h_x_pt = new TH1D("h_x_pt", "h_x_pt", nb_bins_x_tt, xmin_tt, xmax_tt);
    TH1D* h_x_pt_1trk = new TH1D("h_x_pt_1trk", "h_x_pt_1trk", nb_bins_x_tt, xmin_tt, xmax_tt);
    TH1D* h_x_pt_1trk_3pts = new TH1D("h_x_pt_1trk_3pts", "h_x_pt_1trk_3pts", nb_bins_x_tt, xmin_tt, xmax_tt);
    TH1D* h_x_pt_1trk_2pts = new TH1D("h_x_pt_1trk_2pts", "h_x_pt_1trk_2pts", nb_bins_x_tt, xmin_tt, xmax_tt);
    TH1D* h_x_pt_1trk_3pts_dlay = new TH1D("h_x_pt_1trk_3pts_dlay", "h_x_pt_1trk_3pts_dlay", nb_bins_x_tt, xmin_tt, xmax_tt);

    TH1D* h_y_pt = new TH1D("h_y_pt", "h_y_pt", nb_bins_y_tt, ymin_tt, ymax_tt);
    TH1D* h_y_pt_1trk = new TH1D("h_y_pt_1trk", "h_y_pt_1trk", nb_bins_y_tt, ymin_tt, ymax_tt);
    TH1D* h_y_pt_1trk_3pts = new TH1D("h_y_pt_1trk_3pts", "h_y_pt_1trk_3pts", nb_bins_y_tt, ymin_tt, ymax_tt);
    TH1D* h_y_pt_1trk_2pts = new TH1D("h_y_pt_1trk_2pts", "h_y_pt_1trk_2pts", nb_bins_y_tt, ymin_tt, ymax_tt);
    TH1D* h_y_pt_1trk_3pts_dlay = new TH1D("h_y_pt_1trk_3pts_dlay", "h_y_pt_1trk_3pts_dlay", nb_bins_y_tt, ymin_tt, ymax_tt);

    TH1D* h_z_pt = new TH1D("h_z_pt", "h_z_pt", 100, 24.0, 31.0);
    TH1D* h_z_pt_1trk = new TH1D("h_z_pt_1trk", "h_z_pt_1trk", 100, 24.0, 31.0);
    TH1D* h_z_pt_1trk_3pts = new TH1D("h_z_pt_1trk_3pts", "h_z_pt_1trk_3pts", 100, 24.0, 31.0);
    TH1D* h_z_pt_1trk_2pts = new TH1D("h_z_pt_1trk_2pts", "h_z_pt_1trk_2pts", 1000, -24.0, 31.0);
    TH1D* h_z_pt_1trk_3pts_dlay = new TH1D("h_z_pt_1trk_3pts_dlay", "h_z_pt_1trk_3pts_dlay", 100, 24.0, 31.0);

    std::vector<TH2D*> xy_pts_per_layer(6, nullptr);
    for (std::size_t k = 0ul; k < 6ul; ++k) {
        xy_pts_per_layer[k] = new TH2D(Form("xy_pts_per_layer_%zu", k), Form("xy_pts_per_layer_%zu", k), nb_bins_x_tt, xmin_tt, xmax_tt, nb_bins_y_tt, ymin_tt, ymax_tt);
    }
    std::vector<TH2D*> xy_pts_per_layer_3pts(6, nullptr);
    for (std::size_t k = 0ul; k < 6ul; ++k) {
        xy_pts_per_layer_3pts[k] = new TH2D(Form("xy_pts_per_layer_%zu_3pts", k), Form("xy_pts_per_layer_%zu_3pts", k), nb_bins_x_tt, xmin_tt, xmax_tt, nb_bins_y_tt, ymin_tt, ymax_tt);
    }
    std::vector<TH2D*> xy_pts_per_layer_2pts(6, nullptr);
    for (std::size_t k = 0ul; k < 6ul; ++k) {
        xy_pts_per_layer_2pts[k] = new TH2D(Form("xy_pts_per_layer_%zu_2pts", k), Form("xy_pts_per_layer_%zu_2pts", k), nb_bins_x_tt, xmin_tt, xmax_tt, nb_bins_y_tt, ymin_tt, ymax_tt);
    }
    std::vector<TH2D*> xy_pts_per_layer_cdwp_3pts(6, nullptr);
    for (std::size_t k = 0ul; k < 6ul; ++k) {
        xy_pts_per_layer_cdwp_3pts[k] = new TH2D(Form("xy_pts_per_layer_%zu_cdwp_3pts", k), Form("xy_pts_per_layer_%zu_cdwp_3pts", k), nb_bins_x_tt, xmin_tt, xmax_tt, nb_bins_y_tt, ymin_tt, ymax_tt);
    }
    std::vector<TH2D*> xy_pts_per_layer_cdwp_2pts(6, nullptr);
    for (std::size_t k = 0ul; k < 6ul; ++k) {
        xy_pts_per_layer_cdwp_2pts[k] = new TH2D(Form("xy_pts_per_layer_%zu_cdwp_2pts", k), Form("xy_pts_per_layer_%zu_cdwp_2pts", k), nb_bins_x_tt, xmin_tt, xmax_tt, nb_bins_y_tt, ymin_tt, ymax_tt);
    }
    std::vector<TH2D*> xy_pts_per_layer_cdwp_no_tt(3, nullptr);
    for (std::size_t k = 0ul; k < 3ul; ++k) {
        xy_pts_per_layer_cdwp_no_tt[k] = new TH2D(Form("xy_pts_per_layer_%zu_cdwp_no_tt", k), Form("xy_pts_per_layer_%zu_cdwp_no_tt", k), 2 * nb_bins_x_tt, 2 * xmin_tt, 2 * xmax_tt, 2 * nb_bins_y_tt, 2 * ymin_tt, 2 * ymax_tt);
    }
    std::vector<TH2D*> xy_pts_per_layer_cdwp_no_tt_zoom(3, nullptr);
    for (std::size_t k = 0ul; k < 3ul; ++k) {
        xy_pts_per_layer_cdwp_no_tt_zoom[k] = new TH2D(Form("xy_pts_per_layer_%zu_cdwp_no_tt_zoom", k), Form("xy_pts_per_layer_%zu_cdwp_no_tt_zoom", k), nb_bins_x_tt, xmin_tt, xmax_tt, nb_bins_y_tt, ymin_tt, ymax_tt);
    }
    std::vector<TH2D*> xy_pts_per_layer_cdwp_2pts_extrapolation(3, nullptr);
    for (std::size_t k = 0ul; k < 3ul; ++k) {
        xy_pts_per_layer_cdwp_2pts_extrapolation[k] = new TH2D(Form("xy_pts_per_layer_%zu_cdwp_2pts_extrapolation", k), Form("xy_pts_per_layer_%zu_cdwp_2pts_extrapolation", k), nb_bins_x_tt, xmin_tt, xmax_tt, nb_bins_y_tt, ymin_tt, ymax_tt);
    }

    TH1I* h_missing_layer = new TH1I("h_missing_layer", "Missing layer;Layer;Counts", 6, -0.5, 5.5);

    TH1D* h_cdwp_tt_3pts_distance = new TH1D("h_cdwp_tt_3pts_distance", "h_cdwp_tt_3pts_distance", 100, 0.0, 15.0);
    std::vector<TH1D*> cdwp_tt_3pts_distance_per_layer(6, nullptr);
    for (std::size_t k = 0ul; k < 6ul; ++k) {
        cdwp_tt_3pts_distance_per_layer[k] = new TH1D(Form("cdwp_tt_3pts_distance_per_layer_%zu", k), Form("cdwp_tt_3pts_distance_per_layer_%zu", k), 100, 0.0, 15.0);
    }
    TH1D* h_cdwp_tt_2pts_distance = new TH1D("h_cdwp_tt_2pts_distance", "h_cdwp_tt_2pts_distance", 100, 0.0, 15.0);
    std::vector<TH1D*> cdwp_tt_2pts_distance_per_layer(6, nullptr);
    for (std::size_t k = 0ul; k < 6ul; ++k) {
        cdwp_tt_2pts_distance_per_layer[k] = new TH1D(Form("cdwp_tt_2pts_distance_per_layer_%zu", k), Form("cdwp_tt_2pts_distance_per_layer_%zu", k), 100, 0.0, 15.0);
    }

    TH2D* h_cdwp_tt_3pts_x_vs = new TH2D("h_cdwp_tt_3pts_x_vs", "h_cdwp_tt_3pts_x_vs", nb_bins_x_tt, xmin_tt, xmax_tt, nb_bins_x_tt, xmin_tt, xmax_tt);
    TH2D* h_cdwp_tt_3pts_y_vs = new TH2D("h_cdwp_tt_3pts_y_vs", "h_cdwp_tt_3pts_y_vs", nb_bins_y_tt, ymin_tt, ymax_tt, nb_bins_y_tt, ymin_tt, ymax_tt);
    TH2D* h_cdwp_tt_2pts_x_vs = new TH2D("h_cdwp_tt_2pts_x_vs", "h_cdwp_tt_2pts_x_vs", nb_bins_x_tt, xmin_tt, xmax_tt, nb_bins_x_tt, xmin_tt, xmax_tt);
    TH2D* h_cdwp_tt_2pts_y_vs = new TH2D("h_cdwp_tt_2pts_y_vs", "h_cdwp_tt_2pts_y_vs", nb_bins_y_tt, ymin_tt, ymax_tt, nb_bins_y_tt, ymin_tt, ymax_tt);

    TH1D* h_cdwp_tt_3pts_x_diff = new TH1D("h_cdwp_tt_3pts_x_diff", "h_cdwp_tt_3pts_x_diff", 100, -5.0, 5.0);
    std::vector<TH1D*> cdwp_tt_3pts_x_diff_per_layer(6, nullptr);
    for (std::size_t k = 0ul; k < 6ul; ++k) {
        cdwp_tt_3pts_x_diff_per_layer[k] = new TH1D(Form("cdwp_tt_3pts_x_diff_per_layer_%zu", k), Form("cdwp_tt_3pts_x_diff_per_layer_%zu", k), 100, -5.0, 5.0);
    }
    TH1D* h_cdwp_tt_3pts_y_diff = new TH1D("h_cdwp_tt_3pts_y_diff", "h_cdwp_tt_3pts_y_diff", 100, -5.0, 5.0);
    std::vector<TH1D*> cdwp_tt_3pts_y_diff_per_layer(6, nullptr);
    for (std::size_t k = 0ul; k < 6ul; ++k) {
        cdwp_tt_3pts_y_diff_per_layer[k] = new TH1D(Form("cdwp_tt_3pts_y_diff_per_layer_%zu", k), Form("cdwp_tt_3pts_y_diff_per_layer_%zu", k), 100, -5.0, 5.0);
    }
    TH1D* h_cdwp_tt_2pts_x_diff = new TH1D("h_cdwp_tt_2pts_x_diff", "h_cdwp_tt_2pts_x_diff", 100, -5.0, 5.0);
    std::vector<TH1D*> cdwp_tt_2pts_x_diff_per_layer(6, nullptr);
    for (std::size_t k = 0ul; k < 6ul; ++k) {
        cdwp_tt_2pts_x_diff_per_layer[k] = new TH1D(Form("cdwp_tt_2pts_x_diff_per_layer_%zu", k), Form("cdwp_tt_2pts_x_diff_per_layer_%zu", k), 100, -5.0, 5.0);
    }
    TH1D* h_cdwp_tt_2pts_y_diff = new TH1D("h_cdwp_tt_2pts_y_diff", "h_cdwp_tt_2pts_y_diff", 100, -5.0, 5.0);
    std::vector<TH1D*> cdwp_tt_2pts_y_diff_per_layer(6, nullptr);
    for (std::size_t k = 0ul; k < 6ul; ++k) {
        cdwp_tt_2pts_y_diff_per_layer[k] = new TH1D(Form("cdwp_tt_2pts_y_diff_per_layer_%zu", k), Form("cdwp_tt_2pts_y_diff_per_layer_%zu", k), 100, -5.0, 5.0);
    }

    auto layer_id = [](double z) {
        if (24000 <= z && z <= 25000) return 0;  // main
        if (25500 <= z && z <= 26500) return 1;  // main
        if (27000 <= z && z <= 28000) return 2;  // main

        if (30000 <= z && z <= 30200) return 3;  // chimney
        if (30200 <= z && z <= 30400) return 4;  // chimney
        if (30400 <= z && z <= 30600) return 5;  // chimney

        return -1; // not inside any valid layer
    };

    for (int k = 0; k < tree->GetEntries(); ++k) {
        if (k % 1000 == 0) std::cout << "\rEntries: " << k << " / " << tree->GetEntries();
        tree->GetEntry(k);

        if (run_with_tt.find(run_id) == run_with_tt.end()) continue;

        TTimeStamp ts{sec, nsec};

        std::unordered_map<std::string, std::size_t> method_count = {
            {"CdWpTtChi2", 0ul},
            {"WpClassify", 0ul}
        };
        std::size_t j = 0ul;
        for (std::size_t i = 0ul; i < method->size(); ++i) {
            ++(method_count[(*method)[i]]);
            if ((*method)[i] == "CdWpTtChi2") j = i;
        }
        if (/* method_count["WpClassify"] != 1ul || */ method_count["CdWpTtChi2"] != 1ul) continue;

        ipos.SetXYZ((*iposx)[j], (*iposy)[j], (*iposz)[j]);
        fpos.SetXYZ((*fposx)[j], (*fposy)[j], (*fposz)[j]);

        pos_cdwp = ipos;
        dir_cdwp = (fpos - ipos).Unit();

        pos_tt.SetXYZ(Coeff0[0], Coeff1[0], Coeff2[0] + 26452.0);
        dir_tt.SetXYZ(Coeff3[0], Coeff4[0], Coeff5[0]);
        dir_tt = dir_tt.Unit();

        // if (dir_tt.Cross(pos_tt).Mag() / 1000.0 > 17.7) continue;
        
        if (NTotPoints == 0) {
            double layer_tmp[3] = {24500.0, 26000.0, 27500.0};
            for (int i = 0; i < 3; ++i) {
                double t_cdwp = (layer_tmp[i] - pos_cdwp.Z()) / dir_cdwp.Z();
                TVector3 p_cdwp_at_z = pos_cdwp + t_cdwp * dir_cdwp;
                xy_pts_per_layer_cdwp_no_tt[i]->Fill(p_cdwp_at_z.X() / 1000.0, p_cdwp_at_z.Y() / 1000.0);
                xy_pts_per_layer_cdwp_no_tt_zoom[i]->Fill(p_cdwp_at_z.X() / 1000.0, p_cdwp_at_z.Y() / 1000.0);
            }
            continue;
        }

        double ts_diff = 1.0e9 * (ts - *start_TS);

        h_ts_diff->Fill(ts_diff);
        for (int i = 0; i < NTotPoints; ++i) {
            h_x_pt->Fill(PointX[i] / 1000.0);
            h_y_pt->Fill(PointY[i] / 1000.0);
            h_z_pt->Fill((PointZ[i] + 26452.0) / 1000.0);
            int lid = layer_id(PointZ[i] + 26452.0);
            if (lid < 0) continue; 
            xy_pts_per_layer[lid]->Fill(PointX[i] / 1000.0, PointY[i] / 1000.0);
        }
        
        if (NTracks == 1) {

            h_ts_diff_1trk->Fill(ts_diff);
            for (int i = 0; i < NTotPoints; ++i) {
                h_x_pt_1trk->Fill(PointX[i] / 1000.0);
                h_y_pt_1trk->Fill(PointY[i] / 1000.0);
                h_z_pt_1trk->Fill((PointZ[i] + 26452.0) / 1000.0);
            }

            if (NPoints[0] >= 3) {

                h_ts_diff_1trk_3pts->Fill(ts_diff);
                for (int i = 0; i < NTotPoints; ++i) {
                    h_x_pt_1trk_3pts->Fill(PointX[i] / 1000.0);
                    h_y_pt_1trk_3pts->Fill(PointY[i] / 1000.0);
                    h_z_pt_1trk_3pts->Fill((PointZ[i] + 26452.0) / 1000.0);
                    int lid = layer_id(PointZ[i] + 26452.0);
                    if (lid < 0) continue;
                    xy_pts_per_layer_3pts[lid]->Fill(PointX[i] / 1000.0, PointY[i] / 1000.0);
                    double t_cdwp = (PointZ[i] + 26452.0 - pos_cdwp.Z()) / dir_cdwp.Z();
                    TVector3 p_cdwp_at_z = pos_cdwp + t_cdwp * dir_cdwp;
                    double t_tt = (PointZ[i] + 26452.0 - pos_tt.Z()) / dir_tt.Z();
                    TVector3 p_tt_at_z = pos_tt + t_tt * dir_tt;
                    xy_pts_per_layer_cdwp_3pts[lid]->Fill(p_cdwp_at_z.X() / 1000.0, p_cdwp_at_z.Y() / 1000.0);
                    double dist_m = (p_cdwp_at_z - p_tt_at_z).Mag() / 1000.0;
                    h_cdwp_tt_3pts_distance->Fill(dist_m);
                    cdwp_tt_3pts_distance_per_layer[lid]->Fill(dist_m);
                    h_cdwp_tt_3pts_x_vs->Fill(p_cdwp_at_z.X() / 1000.0, PointX[i] / 1000.0);
                    h_cdwp_tt_3pts_y_vs->Fill(p_cdwp_at_z.Y() / 1000.0, PointY[i] / 1000.0);
                    h_cdwp_tt_3pts_x_diff->Fill((p_cdwp_at_z.X() - PointX[i]) / 1000.0);
                    cdwp_tt_3pts_x_diff_per_layer[lid]->Fill((p_cdwp_at_z.X() - PointX[i]) / 1000.0);
                    h_cdwp_tt_3pts_y_diff->Fill((p_cdwp_at_z.Y() - PointY[i]) / 1000.0);
                    cdwp_tt_3pts_y_diff_per_layer[lid]->Fill((p_cdwp_at_z.Y() - PointY[i]) / 1000.0);
                }

                std::unordered_set<int> layers_hit;
                layers_hit.reserve(6);
                for (int i = 0; i < NTotPoints; ++i) {
                    h_z_pt->Fill((PointZ[i] + 26452.0) / 1000.0);
                    int lid = layer_id(PointZ[i] + 26452.0);
                    if (lid < 0) continue;
                    layers_hit.insert(lid);
                }

                if (layers_hit.size() >= 3) { // at least 3 different layers
                    
                    h_ts_diff_1trk_3pts_dlay->Fill(ts_diff);
                    for (int i = 0; i < NTotPoints; ++i) {
                        h_x_pt_1trk_3pts_dlay->Fill(PointX[i] / 1000.0);
                        h_y_pt_1trk_3pts_dlay->Fill(PointY[i] / 1000.0);
                        h_z_pt_1trk_3pts_dlay->Fill((PointZ[i] + 26452.0) / 1000.0);
                    }
                
                }

            }

            else {

                h_ts_diff_1trk_2pts->Fill(ts_diff);
                for (int i = 0; i < NTotPoints; ++i) {
                    h_x_pt_1trk_2pts->Fill(PointX[i] / 1000.0);
                    h_y_pt_1trk_2pts->Fill(PointY[i] / 1000.0);
                    h_z_pt_1trk_2pts->Fill((PointZ[i] + 26452.0) / 1000.0);
                    int lid = layer_id(PointZ[i] + 26452.0);
                    if (lid < 0) continue;
                    xy_pts_per_layer_2pts[lid]->Fill(PointX[i] / 1000.0, PointY[i] / 1000.0);
                    double t_cdwp = (PointZ[i] + 26452.0 - pos_cdwp.Z()) / dir_cdwp.Z();
                    TVector3 p_cdwp_at_z = pos_cdwp + t_cdwp * dir_cdwp;
                    double t_tt = (PointZ[i] + 26452.0 - pos_tt.Z()) / dir_tt.Z();
                    TVector3 p_tt_at_z = pos_tt + t_tt * dir_tt;
                    xy_pts_per_layer_cdwp_2pts[lid]->Fill(p_cdwp_at_z.X() / 1000.0, p_cdwp_at_z.Y() / 1000.0);
                    double dist_m = (p_cdwp_at_z - p_tt_at_z).Mag() / 1000.0;
                    h_cdwp_tt_2pts_distance->Fill(dist_m);
                    cdwp_tt_2pts_distance_per_layer[lid]->Fill(dist_m);
                    h_cdwp_tt_2pts_x_vs->Fill(p_cdwp_at_z.X() / 1000.0, PointX[i] / 1000.0);
                    h_cdwp_tt_2pts_y_vs->Fill(p_cdwp_at_z.Y() / 1000.0, PointY[i] / 1000.0);
                    h_cdwp_tt_2pts_x_diff->Fill((p_cdwp_at_z.X() - PointX[i]) / 1000.0);
                    cdwp_tt_2pts_x_diff_per_layer[lid]->Fill((p_cdwp_at_z.X() - PointX[i]) / 1000.0);
                    h_cdwp_tt_2pts_y_diff->Fill((p_cdwp_at_z.Y() - PointY[i]) / 1000.0);
                    cdwp_tt_2pts_y_diff_per_layer[lid]->Fill((p_cdwp_at_z.Y() - PointY[i]) / 1000.0);
                }

                double layer_tmp[3] = {24500.0, 26000.0, 27500.0};
                for (int i = 0; i < 3; ++i) {
                    double t_cdwp = (layer_tmp[i] - pos_cdwp.Z()) / dir_cdwp.Z();
                    TVector3 p_cdwp_at_z = pos_cdwp + t_cdwp * dir_cdwp;
                    xy_pts_per_layer_cdwp_2pts_extrapolation[i]->Fill(p_cdwp_at_z.X() / 1000.0, p_cdwp_at_z.Y() / 1000.0);
                }

                std::unordered_set<int> layers_hit;
                layers_hit.reserve(6);
                for (int i = 0; i < NTotPoints; ++i) {
                    h_z_pt->Fill((PointZ[i] + 26452.0) / 1000.0);
                    int lid = layer_id(PointZ[i] + 26452.0);
                    if (lid < 0) continue;
                    layers_hit.insert(lid);
                }

                auto find_missings = [&](std::initializer_list<int> group) {
                    int count = 0;
                    int missing = -1;
                    for (int layer : group) {
                        if (layers_hit.count(layer)) count++;
                        else missing = layer;
                    }
                    if (count == group.size() - 1) {
                        h_missing_layer->Fill(missing);
                    }
                };

                find_missings({0, 1, 2});
                find_missings({3, 4, 5});

            }

        }
    }

    auto style = [](TH1* h, Color_t col, int lstyle, int width) {
        h->SetLineColor(col);
        h->SetLineStyle(lstyle);
        h->SetLineWidth(width);
    };

    // ============================================================================================
    // TS diff
    // ============================================================================================

    TCanvas* c_ts_diff = new TCanvas("c_ts_diff", "c_ts_diff", 1000, 1000);
    c_ts_diff->cd();

    style(h_ts_diff, kBlack, kSolid, 3);
    style(h_ts_diff_1trk, kBlue+1, kSolid, 3);
    style(h_ts_diff_1trk_3pts, kRed+1, kSolid, 3);
    style(h_ts_diff_1trk_3pts_dlay, kGreen+2, kSolid, 3);

    h_ts_diff->SetStats(0);
    h_ts_diff->SetMinimum(0.1);
    h_ts_diff->GetXaxis()->SetTitle("TS_{CD} - TS_{TT} (ns)");
    h_ts_diff->GetXaxis()->CenterTitle();
    h_ts_diff->GetYaxis()->SetTitle("Entries");
    h_ts_diff->GetYaxis()->CenterTitle();
    h_ts_diff->GetYaxis()->SetTitleOffset(1.5);
    h_ts_diff->Draw("HIST");
    h_ts_diff_1trk->Draw("HIST SAME");
    h_ts_diff_1trk_3pts->Draw("HIST SAME");
    h_ts_diff_1trk_3pts_dlay->Draw("HIST SAME");

    TLegend* leg_ts_diff = new TLegend(0.15, 0.70, 0.50, 0.85);
    leg_ts_diff->AddEntry(h_ts_diff, "CD/WP + TT correlation", "l");
    leg_ts_diff->AddEntry(h_ts_diff_1trk, "Only 1 track", "l");
    leg_ts_diff->AddEntry(h_ts_diff_1trk_3pts, "At least 3 points", "l");
    leg_ts_diff->AddEntry(h_ts_diff_1trk_3pts_dlay, "Different layers", "l");
    leg_ts_diff->Draw();

    c_ts_diff->SetTickx();
    c_ts_diff->SetTicky();
    c_ts_diff->SetLogy();
    c_ts_diff->Update();

    // ============================================================================================
    // TS fit
    // ============================================================================================

    TCanvas* c_ts_fit = new TCanvas("c_ts_fit", "c_ts_fit", 1000, 1000);
    c_ts_fit->cd();

    TF1* f_model = new TF1("f_model", "[0] + gaus(1)", -1000, 1000);

    f_model->SetParameter(0, 1.0);     // constant
    f_model->SetParameter(1, 500.0);     // gaussian amplitude (rough guess)
    f_model->SetParameter(2, 300.0);     // gaussian mean
    f_model->SetParameter(3, 30.0);      // gaussian sigma

    h_ts_diff->Fit(f_model, "R");

    style(h_ts_diff, kBlack, kSolid, 3);

    h_ts_diff->SetStats(0);
    h_ts_diff->SetMinimum(0.1);
    h_ts_diff->GetXaxis()->SetTitle("TS_{CD} - TS_{TT} (ns)");
    h_ts_diff->GetXaxis()->CenterTitle();
    h_ts_diff->GetYaxis()->SetTitle("Entries");
    h_ts_diff->GetYaxis()->CenterTitle();
    h_ts_diff->GetYaxis()->SetTitleOffset(1.5);
    h_ts_diff->Draw("HIST");

    TF1* f_const = new TF1("f_const", "[0]", -1000, 1000);
    f_const->SetParameter(0, f_model->GetParameter(0));
    f_const->SetLineColor(kBlue);
    f_const->SetLineStyle(kDashed);
    f_const->SetLineWidth(3);

    TF1* f_gaus = new TF1("f_gaus", "gaus(0)", -1000, 1000);
    f_gaus->SetParameter(0, f_model->GetParameter(1)); // amplitude
    f_gaus->SetParameter(1, f_model->GetParameter(2)); // mean
    f_gaus->SetParameter(2, f_model->GetParameter(3)); // sigma
    f_gaus->SetLineColor(kRed);
    f_gaus->SetLineStyle(kDashed);
    f_gaus->SetLineWidth(3);

    f_const->Draw("SAME");
    f_gaus->Draw("SAME");
    f_model->SetLineColor(kRed);
    f_model->SetLineWidth(4);
    f_model->Draw("SAME");

    TLegend* leg_ts_fit = new TLegend(0.15, 0.70, 0.50, 0.85);
    leg_ts_fit->AddEntry(h_ts_diff, "CD/WP + TT correlation", "l");
    leg_ts_fit->AddEntry(f_const, "Constant", "l");
    leg_ts_fit->AddEntry(f_gaus, "Gaussian", "l");
    leg_ts_fit->AddEntry(f_model, "Constant + Gaussian fit", "l");
    leg_ts_fit->Draw();

    c_ts_fit->SetTickx();
    c_ts_fit->SetTicky();
    c_ts_fit->SetLogy();
    c_ts_fit->Update();

    // ============================================================================================
    // TS fit 3 points
    // ============================================================================================

    TCanvas* c_ts_fit_3pts = new TCanvas("c_ts_fit_3pts", "c_ts_fit_3pts", 1000, 1000);
    c_ts_fit_3pts->cd();

    TF1* f_model_3pts = new TF1("f_model_3pts", "gaus(0)", -1000, 1000);

    f_model_3pts->SetParameter(0, 500.0);     // gaussian amplitude (rough guess)
    f_model_3pts->SetParameter(1, 300.0);     // gaussian mean
    f_model_3pts->SetParameter(2, 30.0);      // gaussian sigma

    h_ts_diff_1trk_3pts->Fit(f_model_3pts, "R");

    style(h_ts_diff_1trk_3pts, kBlack, kSolid, 3);

    h_ts_diff_1trk_3pts->SetStats(0);
    h_ts_diff_1trk_3pts->SetMinimum(0.1);
    h_ts_diff_1trk_3pts->GetXaxis()->SetTitle("TS_{CD} - TS_{TT} (ns)");
    h_ts_diff_1trk_3pts->GetXaxis()->CenterTitle();
    h_ts_diff_1trk_3pts->GetYaxis()->SetTitle("Entries");
    h_ts_diff_1trk_3pts->GetYaxis()->CenterTitle();
    h_ts_diff_1trk_3pts->GetYaxis()->SetTitleOffset(1.5);
    h_ts_diff_1trk_3pts->Draw("HIST");

    f_model_3pts->SetLineColor(kRed);
    f_model_3pts->SetLineWidth(4);
    f_model_3pts->Draw("SAME");

    TLegend* leg_ts_fit_3pts = new TLegend(0.15, 0.75, 0.50, 0.85);
    leg_ts_fit_3pts->AddEntry(h_ts_diff, "At least 3 points", "l");
    leg_ts_fit_3pts->AddEntry(f_model_3pts, "Gaussian fit", "l");
    leg_ts_fit_3pts->Draw();

    c_ts_fit_3pts->SetTickx();
    c_ts_fit_3pts->SetTicky();
    c_ts_fit_3pts->SetLogy();
    c_ts_fit_3pts->Update();

    // ============================================================================================
    // TS fit 2 points
    // ============================================================================================

    TCanvas* c_ts_no_corr = new TCanvas("c_ts_no_corr", "c_ts_no_corr", 1000, 1000);
    c_ts_no_corr->cd();

    TF1* f_model_no_corr = new TF1("f_model_no_corr", "[0] + gaus(1)", -1000, 1000);

    f_model_no_corr->SetParameter(0, 1.0);     // constant
    f_model_no_corr->SetParameter(1, 500.0);     // gaussian amplitude (rough guess)
    f_model_no_corr->SetParameter(2, 300.0);     // gaussian mean
    f_model_no_corr->SetParameter(3, 30.0);      // gaussian sigma

    h_ts_diff_1trk_2pts->Fit(f_model_no_corr, "R");

    style(h_ts_diff_1trk_2pts, kBlack, kSolid, 3);

    h_ts_diff_1trk_2pts->SetStats(0);
    h_ts_diff_1trk_2pts->GetXaxis()->SetTitle("TS_{CD} - TS_{TT} (ns)");
    h_ts_diff_1trk_2pts->GetXaxis()->CenterTitle();
    h_ts_diff_1trk_2pts->GetYaxis()->SetTitle("Entries");
    h_ts_diff_1trk_2pts->GetYaxis()->CenterTitle();
    h_ts_diff_1trk_2pts->GetYaxis()->SetTitleOffset(1.5);
    h_ts_diff_1trk_2pts->Draw("HIST");

    TF1* f_const_no_corr = new TF1("f_const_no_corr", "[0]", -1000, 1000);
    f_const_no_corr->SetParameter(0, f_model_no_corr->GetParameter(0));
    f_const_no_corr->SetLineColor(kBlue);
    f_const_no_corr->SetLineStyle(kDashed);
    f_const_no_corr->SetLineWidth(3);

    TF1* f_gaus_no_corr = new TF1("f_gaus_no_corr", "gaus(0)", -1000, 1000);
    f_gaus_no_corr->SetParameter(0, f_model_no_corr->GetParameter(1)); // amplitude
    f_gaus_no_corr->SetParameter(1, f_model_no_corr->GetParameter(2)); // mean
    f_gaus_no_corr->SetParameter(2, f_model_no_corr->GetParameter(3)); // sigma
    f_gaus_no_corr->SetLineColor(kRed);
    f_gaus_no_corr->SetLineStyle(kDashed);
    f_gaus_no_corr->SetLineWidth(3);

    f_const_no_corr->Draw("SAME");
    f_gaus_no_corr->Draw("SAME");
    f_model_no_corr->SetLineColor(kRed);
    f_model_no_corr->SetLineWidth(4);
    f_model_no_corr->Draw("SAME");

    TLegend* leg_ts_no_corr = new TLegend(0.15, 0.70, 0.50, 0.85);
    leg_ts_no_corr->AddEntry(h_ts_diff_1trk_2pts, "At most 2 points", "l");
    leg_ts_no_corr->AddEntry(f_const_no_corr, "Constant", "l");
    leg_ts_no_corr->AddEntry(f_gaus_no_corr, "Gaussian", "l");
    leg_ts_no_corr->AddEntry(f_model_no_corr, "Constant + Gaussian fit", "l");
    leg_ts_no_corr->Draw();

    c_ts_no_corr->SetTickx();
    c_ts_no_corr->SetTicky();
    c_ts_no_corr->SetLogy();
    c_ts_no_corr->Update();

    // ============================================================================================
    // TS 3 points vs 2 points
    // ============================================================================================

    TCanvas* c_ts_3_vs_2 = new TCanvas("c_ts_3_vs_2", "c_ts_3_vs_2", 1000, 1000);
    c_ts_3_vs_2->cd();

    style(h_ts_diff_1trk_3pts, kRed+1, kSolid, 3);
    style(h_ts_diff_1trk_2pts, kViolet-3, kDashed, 3);

    h_ts_diff_1trk_3pts->SetStats(0);
    h_ts_diff_1trk_3pts->GetXaxis()->SetTitle("TS_{CD} - TS_{TT} (ns)");
    h_ts_diff_1trk_3pts->GetXaxis()->CenterTitle();
    h_ts_diff_1trk_3pts->GetYaxis()->SetTitle("Entries");
    h_ts_diff_1trk_3pts->GetYaxis()->CenterTitle();
    h_ts_diff_1trk_3pts->GetYaxis()->SetTitleOffset(1.5);
    h_ts_diff_1trk_3pts->Draw("HIST");
    h_ts_diff_1trk_2pts->Draw("HIST SAME");

    std::cout << "Number of 1 track 3 points = " << h_ts_diff_1trk_3pts->GetEntries() << '\n';
    std::cout << "Number of 1 track 2 points = " << h_ts_diff_1trk_2pts->GetEntries() << '\n';

    TLegend* leg_ts_3vs_2 = new TLegend(0.15, 0.75, 0.50, 0.85);
    leg_ts_3vs_2->AddEntry(h_ts_diff_1trk_3pts, "At least 3 points", "l");
    leg_ts_3vs_2->AddEntry(h_ts_diff_1trk_2pts, "At most 2 points", "l");
    leg_ts_3vs_2->Draw();

    c_ts_3_vs_2->SetTickx();
    c_ts_3_vs_2->SetTicky();
    c_ts_3_vs_2->SetLogy();
    c_ts_3_vs_2->Update();

    // ============================================================================================
    // Point X - 3 points
    // ============================================================================================

    TCanvas* c_x_pt = new TCanvas("c_x_pt", "c_x_pt", 1000, 1000);
    c_x_pt->cd();

    style(h_x_pt, kBlack, kSolid, 3);
    style(h_x_pt_1trk, kBlue+1, kSolid, 3);
    style(h_x_pt_1trk_3pts, kRed+1, kSolid, 3);
    style(h_x_pt_1trk_3pts_dlay, kGreen+2, kSolid, 3);

    h_x_pt->SetStats(0);
    h_x_pt->SetMinimum(0.1);
    h_x_pt->GetXaxis()->SetTitle("x (m)");
    h_x_pt->GetXaxis()->CenterTitle(kTRUE);
    h_x_pt->GetYaxis()->SetTitle("Entries");
    h_x_pt->GetYaxis()->CenterTitle(kTRUE);
    h_x_pt->GetYaxis()->SetTitleOffset(1.25);
    h_x_pt->Draw("HIST");
    h_x_pt_1trk->Draw("HIST SAME");
    h_x_pt_1trk_3pts->Draw("HIST SAME");
    h_x_pt_1trk_3pts_dlay->Draw("HIST SAME");

    TLegend* leg_x_pt = new TLegend(0.35, 0.15, 0.65, 0.30);
    leg_x_pt->AddEntry(h_x_pt, "CD/WP + TT correlation", "l");
    leg_x_pt->AddEntry(h_x_pt_1trk, "Only 1 track", "l");
    leg_x_pt->AddEntry(h_x_pt_1trk_3pts, "At least 3 points", "l");
    leg_x_pt->AddEntry(h_x_pt_1trk_3pts_dlay, "Different layers", "l");
    leg_x_pt->Draw();

    c_x_pt->SetTickx();
    c_x_pt->SetTicky();
    c_x_pt->SetLogy();
    c_x_pt->Update();

    // ============================================================================================
    // Point Y - 3 points
    // ============================================================================================

    TCanvas* c_y_pt = new TCanvas("c_y_pt", "c_y_pt", 1000, 1000);
    c_y_pt->cd();

    style(h_y_pt, kBlack, kSolid, 3);
    style(h_y_pt_1trk, kBlue+1, kSolid, 3);
    style(h_y_pt_1trk_3pts, kRed+1, kSolid, 3);
    style(h_y_pt_1trk_3pts_dlay, kGreen+2, kSolid, 3);

    h_y_pt->SetStats(0);
    h_y_pt->SetMinimum(0.1);
    h_y_pt->GetXaxis()->SetTitle("y (m)");
    h_y_pt->GetXaxis()->CenterTitle(kTRUE);
    h_y_pt->GetYaxis()->SetTitle("Entries");
    h_y_pt->GetYaxis()->CenterTitle(kTRUE);
    h_y_pt->GetYaxis()->SetTitleOffset(1.25);
    h_y_pt->Draw("HIST");
    h_y_pt_1trk->Draw("HIST SAME");
    h_y_pt_1trk_3pts->Draw("HIST SAME");
    h_y_pt_1trk_3pts_dlay->Draw("HIST SAME");

    TLegend* leg_y_pt = new TLegend(0.35, 0.15, 0.65, 0.30);
    leg_y_pt->AddEntry(h_y_pt, "CD/WP + TT correlation", "l");
    leg_y_pt->AddEntry(h_y_pt_1trk, "Only 1 track", "l");
    leg_y_pt->AddEntry(h_y_pt_1trk_3pts, "At least 3 points", "l");
    leg_y_pt->AddEntry(h_y_pt_1trk_3pts_dlay, "Different layers", "l");
    leg_y_pt->Draw();

    c_y_pt->SetTickx();
    c_y_pt->SetTicky();
    c_y_pt->SetLogy();
    c_y_pt->Update();

    // ============================================================================================
    // Point Z - 3 points
    // ============================================================================================

    TCanvas* c_z_pt = new TCanvas("c_z_pt", "c_z_pt", 1000, 1000);
    c_z_pt->cd();

    style(h_z_pt, kBlack, kSolid, 3);
    style(h_z_pt_1trk, kBlue+1, kSolid, 3);
    style(h_z_pt_1trk_3pts, kRed+1, kSolid, 3);
    style(h_z_pt_1trk_3pts_dlay, kGreen+2, kSolid, 3);

    h_z_pt->SetStats(0);
    h_z_pt->SetMinimum(0.1);
    h_z_pt->GetXaxis()->SetTitle("z (m)");
    h_z_pt->GetXaxis()->CenterTitle(kTRUE);
    h_z_pt->GetYaxis()->SetTitle("Entries");
    h_z_pt->GetYaxis()->CenterTitle(kTRUE);
    h_z_pt->GetYaxis()->SetTitleOffset(1.25);
    h_z_pt->Draw("HIST");
    h_z_pt_1trk->Draw("HIST SAME");
    h_z_pt_1trk_3pts->Draw("HIST SAME");
    h_z_pt_1trk_3pts_dlay->Draw("HIST SAME");

    TLegend* leg_z_pt = new TLegend(0.60, 0.70, 0.85, 0.85);
    leg_z_pt->AddEntry(h_z_pt, "CD/WP + TT correlation", "l");
    leg_z_pt->AddEntry(h_z_pt_1trk, "Only 1 track", "l");
    leg_z_pt->AddEntry(h_z_pt_1trk_3pts, "At least 3 points", "l");
    leg_z_pt->AddEntry(h_z_pt_1trk_3pts_dlay, "Different layers", "l");
    leg_z_pt->Draw();

    c_z_pt->SetTickx();
    c_z_pt->SetTicky();
    c_z_pt->SetLogy();
    c_z_pt->Update();

    // ============================================================================================
    // Point X - 3 vs 2 points
    // ============================================================================================

    TCanvas* c_x_pt_3_vs_2 = new TCanvas("c_x_pt_3_vs_2", "c_x_pt_3_vs_2", 1000, 1000);
    c_x_pt_3_vs_2->cd();

    style(h_x_pt_1trk_3pts, kRed+1, kSolid, 3);
    style(h_x_pt_1trk_2pts, kViolet-3, kDashed, 3);

    h_x_pt_1trk_3pts->SetStats(0);
    h_x_pt_1trk_3pts->SetMinimum(0.1);
    h_x_pt_1trk_3pts->GetXaxis()->SetTitle("x (m)");
    h_x_pt_1trk_3pts->GetXaxis()->CenterTitle(kTRUE);
    h_x_pt_1trk_3pts->GetYaxis()->SetTitle("Entries");
    h_x_pt_1trk_3pts->GetYaxis()->CenterTitle(kTRUE);
    h_x_pt_1trk_3pts->GetYaxis()->SetTitleOffset(1.25);
    h_x_pt_1trk_3pts->Draw("HIST");
    h_x_pt_1trk_2pts->Draw("HIST SAME");

    TLegend* leg_x_pt_3_vs_2 = new TLegend(0.35, 0.15, 0.65, 0.25);
    leg_x_pt_3_vs_2->AddEntry(h_x_pt_1trk_3pts, "At least 3 points", "l");
    leg_x_pt_3_vs_2->AddEntry(h_x_pt_1trk_2pts, "At most 2 points", "l");
    leg_x_pt_3_vs_2->Draw();

    c_x_pt_3_vs_2->SetTickx();
    c_x_pt_3_vs_2->SetTicky();
    c_x_pt_3_vs_2->SetLogy();
    c_x_pt_3_vs_2->Update();

    // ============================================================================================
    // Point Y - 3 vs 2 points
    // ============================================================================================

    TCanvas* c_y_pt_3_vs_2 = new TCanvas("c_y_pt_3_vs_2", "c_y_pt_3_vs_2", 1000, 1000);
    c_y_pt_3_vs_2->cd();

    style(h_y_pt_1trk_3pts, kRed+1, kSolid, 3);
    style(h_y_pt_1trk_2pts, kViolet-3, kDashed, 3);

    h_y_pt_1trk_3pts->SetStats(0);
    h_y_pt_1trk_3pts->SetMinimum(0.1);
    h_y_pt_1trk_3pts->GetXaxis()->SetTitle("y (m)");
    h_y_pt_1trk_3pts->GetXaxis()->CenterTitle(kTRUE);
    h_y_pt_1trk_3pts->GetYaxis()->SetTitle("Entries");
    h_y_pt_1trk_3pts->GetYaxis()->CenterTitle(kTRUE);
    h_y_pt_1trk_3pts->GetYaxis()->SetTitleOffset(1.25);
    h_y_pt_1trk_3pts->Draw("HIST");
    h_y_pt_1trk_2pts->Draw("HIST SAME");

    TLegend* leg_y_pt_3_vs_2 = new TLegend(0.35, 0.15, 0.65, 0.25);
    leg_y_pt_3_vs_2->AddEntry(h_y_pt_1trk_3pts, "At least 3 points", "l");
    leg_y_pt_3_vs_2->AddEntry(h_y_pt_1trk_2pts, "At most 2 points", "l");
    leg_y_pt_3_vs_2->Draw();

    c_y_pt_3_vs_2->SetTickx();
    c_y_pt_3_vs_2->SetTicky();
    c_y_pt_3_vs_2->SetLogy();
    c_y_pt_3_vs_2->Update();

    // ============================================================================================
    // Point Z - 3 vs 2 points
    // ============================================================================================

    TCanvas* c_z_pt_3_vs_2 = new TCanvas("c_z_pt_3_vs_2", "c_z_pt_3_vs_2", 1000, 1000);
    c_z_pt_3_vs_2->cd();

    style(h_z_pt_1trk_3pts, kRed+1, kSolid, 3);
    style(h_z_pt_1trk_2pts, kViolet-3, kDashed, 3);

    h_z_pt_1trk_3pts->SetStats(0);
    h_z_pt_1trk_3pts->SetMinimum(0.1);
    h_z_pt_1trk_3pts->GetXaxis()->SetTitle("z (m)");
    h_z_pt_1trk_3pts->GetXaxis()->CenterTitle(kTRUE);
    h_z_pt_1trk_3pts->GetYaxis()->SetTitle("Entries");
    h_z_pt_1trk_3pts->GetYaxis()->CenterTitle(kTRUE);
    h_z_pt_1trk_3pts->GetYaxis()->SetTitleOffset(1.25);
    h_z_pt_1trk_3pts->Draw("HIST");
    h_z_pt_1trk_2pts->Draw("HIST SAME");

    TLegend* leg_z_pt_3_vs_2 = new TLegend(0.60, 0.75, 0.85, 0.85);
    leg_z_pt_3_vs_2->AddEntry(h_z_pt_1trk_3pts, "At least 3 points", "l");
    leg_z_pt_3_vs_2->AddEntry(h_z_pt_1trk_2pts, "At most 2 points", "l");
    leg_z_pt_3_vs_2->Draw();

    c_z_pt_3_vs_2->SetTickx();
    c_z_pt_3_vs_2->SetTicky();
    c_z_pt_3_vs_2->SetLogy();
    c_z_pt_3_vs_2->Update();

    // ============================================================================================
    // X-Y map for each main layer
    // ============================================================================================

    TCanvas* c_xy_stack = new TCanvas("c_xy_stack", "XY Points per Layer", 800, 1000);
    c_xy_stack->Divide(1, 6);

    const int layers[6] = {5, 4, 3, 2, 1, 0};

    for (int p = 0; p < 6; p++) {
        TVirtualPad* pad = c_xy_stack->cd(p + 1);
        pad->SetBottomMargin(0.20);

        TH2D* h = xy_pts_per_layer[layers[p]];
        h->SetStats(0);
        h->SetTitle("");
        h->GetXaxis()->SetLabelSize(0.12);
        h->GetYaxis()->SetLabelSize(0.12);
        h->GetZaxis()->SetLabelSize(0.12);
        h->GetXaxis()->SetTitleSize(0.12);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetZaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetTitleOffset(0.75);
        h->GetYaxis()->SetTitleOffset(0.25);
        h->GetZaxis()->SetTitleOffset(0.15);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitle("x (m)");
        h->GetYaxis()->SetTitle("y (m)");
        h->GetZaxis()->SetTitle("Entries");
        h->Draw("COLZ");

        pad->SetTickx();
        pad->SetTicky();
    }

    c_xy_stack->Update();

    // ============================================================================================
    // X-Y map for each main layer - 3 points
    // ============================================================================================

    TCanvas* c_xy_3pts_stack = new TCanvas("c_xy_3pts_stack", "XY Points per Layer - 3 points", 800, 1000);
    c_xy_3pts_stack->Divide(1, 6);

    for (int p = 0; p < 6; p++) {
        TVirtualPad* pad = c_xy_3pts_stack->cd(p + 1);
        pad->SetBottomMargin(0.20);

        TH2D* h = xy_pts_per_layer_3pts[layers[p]];
        h->SetStats(0);
        h->SetTitle("");
        h->GetXaxis()->SetLabelSize(0.12);
        h->GetYaxis()->SetLabelSize(0.12);
        h->GetZaxis()->SetLabelSize(0.12);
        h->GetXaxis()->SetTitleSize(0.12);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetZaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetTitleOffset(0.75);
        h->GetYaxis()->SetTitleOffset(0.25);
        h->GetZaxis()->SetTitleOffset(0.15);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitle("x (m)");
        h->GetYaxis()->SetTitle("y (m)");
        h->GetZaxis()->SetTitle("Entries");
        h->Draw("COLZ");

        pad->SetTickx();
        pad->SetTicky();
    }

    c_xy_3pts_stack->Update();

    // ============================================================================================
    // X-Y map for each main layer - 2 points
    // ============================================================================================

    TCanvas* c_xy_2pts_stack = new TCanvas("c_xy_2pts_stack", "XY Points per Layer - 2 Points", 800, 1000);
    c_xy_2pts_stack->Divide(1, 6);

    for (int p = 0; p < 6; p++) {
        TVirtualPad* pad = c_xy_2pts_stack->cd(p + 1);
        pad->SetBottomMargin(0.20);

        TH2D* h = xy_pts_per_layer_2pts[layers[p]];
        h->SetStats(0);
        h->SetTitle("");
        h->GetXaxis()->SetLabelSize(0.12);
        h->GetYaxis()->SetLabelSize(0.12);
        h->GetZaxis()->SetLabelSize(0.12);
        h->GetXaxis()->SetTitleSize(0.12);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetZaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetTitleOffset(0.75);
        h->GetYaxis()->SetTitleOffset(0.25);
        h->GetZaxis()->SetTitleOffset(0.15);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitle("x (m)");
        h->GetYaxis()->SetTitle("y (m)");
        h->GetZaxis()->SetTitle("Entries");
        h->Draw("COLZ");

        pad->SetTickx();
        pad->SetTicky();
    }

    c_xy_2pts_stack->Update();

    // ============================================================================================
    // X-Y map for each main layer - 3 points
    // ============================================================================================

    TCanvas* c_xy_cdwp_3pts_stack = new TCanvas("c_xy_cdwp_3pts_stack", "XY Points per Layer - CD/WP 3 points", 800, 1000);
    c_xy_cdwp_3pts_stack->Divide(1, 6);

    for (int p = 0; p < 6; p++) {
        TVirtualPad* pad = c_xy_cdwp_3pts_stack->cd(p + 1);
        pad->SetBottomMargin(0.20);

        TH2D* h = xy_pts_per_layer_cdwp_3pts[layers[p]];
        h->SetStats(0);
        h->SetTitle("");
        h->GetXaxis()->SetLabelSize(0.12);
        h->GetYaxis()->SetLabelSize(0.12);
        h->GetZaxis()->SetLabelSize(0.12);
        h->GetXaxis()->SetTitleSize(0.12);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetZaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetTitleOffset(0.75);
        h->GetYaxis()->SetTitleOffset(0.25);
        h->GetZaxis()->SetTitleOffset(0.15);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitle("x (m)");
        h->GetYaxis()->SetTitle("y (m)");
        h->GetZaxis()->SetTitle("Entries");
        h->Draw("COLZ");

        pad->SetTickx();
        pad->SetTicky();
    }

    c_xy_cdwp_3pts_stack->Update();

    // ============================================================================================
    // X-Y map for each main layer - 2 points
    // ============================================================================================

    TCanvas* c_xy_cdwp_2pts_stack = new TCanvas("c_xy_cdwp_2pts_stack", "XY Points per Layer - CD/WP 2 points", 800, 1000);
    c_xy_cdwp_2pts_stack->Divide(1, 6);

    for (int p = 0; p < 6; p++) {
        TVirtualPad* pad = c_xy_cdwp_2pts_stack->cd(p + 1);
        pad->SetBottomMargin(0.20);

        TH2D* h = xy_pts_per_layer_cdwp_2pts[layers[p]];
        h->SetStats(0);
        h->SetTitle("");
        h->GetXaxis()->SetLabelSize(0.12);
        h->GetYaxis()->SetLabelSize(0.12);
        h->GetZaxis()->SetLabelSize(0.12);
        h->GetXaxis()->SetTitleSize(0.12);
        h->GetYaxis()->SetTitleSize(0.12);
        h->GetZaxis()->SetTitleSize(0.12);
        h->GetXaxis()->SetTitleOffset(0.75);
        h->GetYaxis()->SetTitleOffset(0.25);
        h->GetZaxis()->SetTitleOffset(0.15);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitle("x (m)");
        h->GetYaxis()->SetTitle("y (m)");
        h->GetZaxis()->SetTitle("Entries");
        h->Draw("COLZ");

        pad->SetTickx();
        pad->SetTicky();
    }

    c_xy_cdwp_2pts_stack->Update();

    // ============================================================================================
    // Distance between CD/WP and TT points - 3 points
    // ============================================================================================

    TCanvas* c_dist_cdwp_tt_3pts = new TCanvas("c_dist_cdwp_tt_3pts", "c_dist_cdwp_tt_3pts", 1000, 1000);
    c_dist_cdwp_tt_3pts->cd();

    style(h_cdwp_tt_3pts_distance, kBlack, kSolid, 3);
    style(cdwp_tt_3pts_distance_per_layer[0], kBlue+1, kDashed, 3);
    style(cdwp_tt_3pts_distance_per_layer[1], kRed+1, kDashed, 3);
    style(cdwp_tt_3pts_distance_per_layer[2], kGreen+2, kDashed, 3);

    h_cdwp_tt_3pts_distance->SetMinimum(0.1);
    h_cdwp_tt_3pts_distance->GetXaxis()->SetTitle("d (m)");
    h_cdwp_tt_3pts_distance->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_3pts_distance->GetYaxis()->SetTitle("Entries");
    h_cdwp_tt_3pts_distance->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_3pts_distance->GetYaxis()->SetTitleOffset(1.25);

    h_cdwp_tt_3pts_distance->Draw("HIST");
    cdwp_tt_3pts_distance_per_layer[0]->Draw("HIST SAME");
    cdwp_tt_3pts_distance_per_layer[1]->Draw("HIST SAME");
    cdwp_tt_3pts_distance_per_layer[2]->Draw("HIST SAME");

    TLegend* leg_cdwp_tt_3pts_distance = new TLegend(0.25, 0.70, 0.45, 0.85);
    leg_cdwp_tt_3pts_distance->AddEntry(h_cdwp_tt_3pts_distance, "All", "l");
    leg_cdwp_tt_3pts_distance->AddEntry(cdwp_tt_3pts_distance_per_layer[0], "Layer 1", "l");
    leg_cdwp_tt_3pts_distance->AddEntry(cdwp_tt_3pts_distance_per_layer[1], "Layer 2", "l");
    leg_cdwp_tt_3pts_distance->AddEntry(cdwp_tt_3pts_distance_per_layer[2], "Layer 3", "l");
    leg_cdwp_tt_3pts_distance->Draw();

    c_dist_cdwp_tt_3pts->SetTickx();
    c_dist_cdwp_tt_3pts->SetTicky();
    c_dist_cdwp_tt_3pts->SetLogy();
    c_dist_cdwp_tt_3pts->Update();

    TPaveStats* st_cdwp_tt_3pts_distance = (TPaveStats*)h_cdwp_tt_3pts_distance->FindObject("stats");
    st_cdwp_tt_3pts_distance->SetOptStat(ToROOTOpt(StatOpt::Entries | StatOpt::Mean | StatOpt::RMS | StatOpt::Underflow | StatOpt::Overflow));
    st_cdwp_tt_3pts_distance->SetX1NDC(0.5);
    st_cdwp_tt_3pts_distance->SetX2NDC(0.85);
    st_cdwp_tt_3pts_distance->SetY1NDC(0.70);
    st_cdwp_tt_3pts_distance->SetY2NDC(0.85);

    c_dist_cdwp_tt_3pts->Modified();
    c_dist_cdwp_tt_3pts->Update();

    // ============================================================================================
    // Distance between CD/WP and TT points - 2 points
    // ============================================================================================

    TCanvas* c_dist_cdwp_tt_2pts = new TCanvas("c_dist_cdwp_tt_2pts", "c_dist_cdwp_tt_2pts", 1000, 1000);
    c_dist_cdwp_tt_2pts->cd();

    style(h_cdwp_tt_2pts_distance, kBlack, kSolid, 3);
    style(cdwp_tt_2pts_distance_per_layer[0], kBlue+1, kDashed, 3);
    style(cdwp_tt_2pts_distance_per_layer[1], kRed+1, kDashed, 3);
    style(cdwp_tt_2pts_distance_per_layer[2], kGreen+2, kDashed, 3);

    h_cdwp_tt_2pts_distance->SetMinimum(0.1);
    h_cdwp_tt_2pts_distance->GetXaxis()->SetTitle("d (m)");
    h_cdwp_tt_2pts_distance->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_2pts_distance->GetYaxis()->SetTitle("Entries");
    h_cdwp_tt_2pts_distance->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_2pts_distance->GetYaxis()->SetTitleOffset(1.25);

    h_cdwp_tt_2pts_distance->Draw("HIST");
    cdwp_tt_2pts_distance_per_layer[0]->Draw("HIST SAME");
    cdwp_tt_2pts_distance_per_layer[1]->Draw("HIST SAME");
    cdwp_tt_2pts_distance_per_layer[2]->Draw("HIST SAME");

    TLegend* leg_cdwp_tt_2pts_distance = new TLegend(0.25, 0.70, 0.45, 0.85);
    leg_cdwp_tt_2pts_distance->AddEntry(h_cdwp_tt_2pts_distance, "All", "l");
    leg_cdwp_tt_2pts_distance->AddEntry(cdwp_tt_2pts_distance_per_layer[0], "Layer 1", "l");
    leg_cdwp_tt_2pts_distance->AddEntry(cdwp_tt_2pts_distance_per_layer[1], "Layer 2", "l");
    leg_cdwp_tt_2pts_distance->AddEntry(cdwp_tt_2pts_distance_per_layer[2], "Layer 3", "l");
    leg_cdwp_tt_2pts_distance->Draw();

    c_dist_cdwp_tt_2pts->SetTickx();
    c_dist_cdwp_tt_2pts->SetTicky();
    c_dist_cdwp_tt_2pts->SetLogy();
    c_dist_cdwp_tt_2pts->Update();

    TPaveStats* st_cdwp_tt_2pts_distance = (TPaveStats*)h_cdwp_tt_2pts_distance->FindObject("stats");
    st_cdwp_tt_2pts_distance->SetOptStat(ToROOTOpt(StatOpt::Entries | StatOpt::Mean | StatOpt::RMS | StatOpt::Underflow | StatOpt::Overflow));
    st_cdwp_tt_2pts_distance->SetX1NDC(0.5);
    st_cdwp_tt_2pts_distance->SetX2NDC(0.85);
    st_cdwp_tt_2pts_distance->SetY1NDC(0.70);
    st_cdwp_tt_2pts_distance->SetY2NDC(0.85);

    c_dist_cdwp_tt_2pts->Modified();
    c_dist_cdwp_tt_2pts->Update();

    // ============================================================================================
    // CD/WP x vs TT x - 3 points
    // ============================================================================================

    TCanvas* c_cdwp_tt_3pts_x_vs = new TCanvas("c_cdwp_tt_3pts_x_vs", "c_cdwp_tt_3pts_x_vs", 1000, 1000);
    c_cdwp_tt_3pts_x_vs->cd();

    h_cdwp_tt_3pts_x_vs->SetStats(0);
    h_cdwp_tt_3pts_x_vs->GetXaxis()->SetTitle("x (m) [CD/WP]");
    h_cdwp_tt_3pts_x_vs->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_3pts_x_vs->GetYaxis()->SetTitle("x (m) [TT]");
    h_cdwp_tt_3pts_x_vs->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_3pts_x_vs->GetYaxis()->SetTitleOffset(1.25);
    h_cdwp_tt_3pts_x_vs->GetZaxis()->SetTitle("Entries");
    h_cdwp_tt_3pts_x_vs->Draw("COLZ");

    c_cdwp_tt_3pts_x_vs->SetTickx();
    c_cdwp_tt_3pts_x_vs->SetTicky();
    c_cdwp_tt_3pts_x_vs->Update();

    // ============================================================================================
    // CD/WP y vs TT y - 3 points
    // ============================================================================================

    TCanvas* c_cdwp_tt_3pts_y_vs = new TCanvas("c_cdwp_tt_3pts_y_vs", "c_cdwp_tt_3pts_y_vs", 1000, 1000);
    c_cdwp_tt_3pts_y_vs->cd();

    h_cdwp_tt_3pts_y_vs->SetStats(0);
    h_cdwp_tt_3pts_y_vs->GetXaxis()->SetTitle("y (m) [CD/WP]");
    h_cdwp_tt_3pts_y_vs->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_3pts_y_vs->GetYaxis()->SetTitle("y (m) [TT]");
    h_cdwp_tt_3pts_y_vs->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_3pts_y_vs->GetYaxis()->SetTitleOffset(1.25);
    h_cdwp_tt_3pts_y_vs->GetZaxis()->SetTitle("Entries");
    h_cdwp_tt_3pts_y_vs->Draw("COLZ");

    c_cdwp_tt_3pts_y_vs->SetTickx();
    c_cdwp_tt_3pts_y_vs->SetTicky();
    c_cdwp_tt_3pts_y_vs->Update();

    // ============================================================================================
    // CD/WP x vs TT x - 2 points
    // ============================================================================================

    TCanvas* c_cdwp_tt_2pts_x_vs = new TCanvas("c_cdwp_tt_2pts_x_vs", "c_cdwp_tt_2pts_x_vs", 1000, 1000);
    c_cdwp_tt_2pts_x_vs->cd();

    h_cdwp_tt_2pts_x_vs->SetStats(0);
    h_cdwp_tt_2pts_x_vs->GetXaxis()->SetTitle("x (m) [CD/WP]");
    h_cdwp_tt_2pts_x_vs->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_2pts_x_vs->GetYaxis()->SetTitle("x (m) [TT]");
    h_cdwp_tt_2pts_x_vs->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_2pts_x_vs->GetYaxis()->SetTitleOffset(1.25);
    h_cdwp_tt_2pts_x_vs->GetZaxis()->SetTitle("Entries");
    h_cdwp_tt_2pts_x_vs->Draw("COLZ");

    c_cdwp_tt_2pts_x_vs->SetTickx();
    c_cdwp_tt_2pts_x_vs->SetTicky();
    c_cdwp_tt_2pts_x_vs->Update();

    // ============================================================================================
    // CD/WP y vs TT y - 2 points
    // ============================================================================================

    TCanvas* c_cdwp_tt_2pts_y_vs = new TCanvas("c_cdwp_tt_2pts_y_vs", "c_cdwp_tt_2pts_y_vs", 1000, 1000);
    c_cdwp_tt_2pts_y_vs->cd();

    h_cdwp_tt_2pts_y_vs->SetStats(0);
    h_cdwp_tt_2pts_y_vs->GetXaxis()->SetTitle("y (m) [CD/WP]");
    h_cdwp_tt_2pts_y_vs->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_2pts_y_vs->GetYaxis()->SetTitle("y (m) [TT]");
    h_cdwp_tt_2pts_y_vs->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_2pts_y_vs->GetYaxis()->SetTitleOffset(1.25);
    h_cdwp_tt_2pts_y_vs->GetZaxis()->SetTitle("Entries");
    h_cdwp_tt_2pts_y_vs->Draw("COLZ");

    c_cdwp_tt_2pts_y_vs->SetTickx();
    c_cdwp_tt_2pts_y_vs->SetTicky();
    c_cdwp_tt_2pts_y_vs->Update();

    // ============================================================================================
    // (CD/WP - TT) - x - 3 points
    // ============================================================================================

    TCanvas* c_cdwp_tt_3pts_x_diff = new TCanvas("c_cdwp_tt_3pts_x_diff", "c_cdwp_tt_3pts_x_diff", 1000, 1000);
    c_cdwp_tt_3pts_x_diff->cd();

    style(h_cdwp_tt_3pts_x_diff, kBlack, kSolid, 3);
    style(cdwp_tt_3pts_x_diff_per_layer[0], kBlue+1, kDashed, 3);
    style(cdwp_tt_3pts_x_diff_per_layer[1], kRed+1, kDashed, 3);
    style(cdwp_tt_3pts_x_diff_per_layer[2], kGreen-2, kDashed, 3);

    // TF1* laplace_cdwp_tt_3pts_x_diff = new TF1("laplace_cdwp_tt_3pts_x_diff","[0] * exp(-abs(x - [1]) / [2]) + [3]", -5.0, 5.0);
    // laplace_cdwp_tt_3pts_x_diff->SetParameters(100.0, 0.0, 0.1, 1.0); // [0]=A, [1]=mu, [2]=tau, [3]=C
    TF1* laplace_cdwp_tt_3pts_x_diff = new TF1("laplace_cdwp_tt_3pts_x_diff","[0]*exp(-abs(x-[1])/[2]) + [3]*exp(-abs(x-[4])/[5])", -5.0, 5.0);
    laplace_cdwp_tt_3pts_x_diff->SetParameters(
        h_cdwp_tt_3pts_x_diff->GetMaximum(), 0.0, 0.5, 
        0.25 * h_cdwp_tt_3pts_x_diff->GetMaximum(), 0.0, 0.25
    );
    
    h_cdwp_tt_3pts_x_diff->Fit(laplace_cdwp_tt_3pts_x_diff, "R");
    h_cdwp_tt_3pts_x_diff->GetXaxis()->SetTitle("#Delta x (m) [CD/WP - TT]");
    h_cdwp_tt_3pts_x_diff->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_3pts_x_diff->GetYaxis()->SetTitle("Entries");
    h_cdwp_tt_3pts_x_diff->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_3pts_x_diff->GetYaxis()->SetTitleOffset(1.25);
    h_cdwp_tt_3pts_x_diff->Draw("HIST");
    cdwp_tt_3pts_x_diff_per_layer[0]->Draw("HIST SAME");
    cdwp_tt_3pts_x_diff_per_layer[1]->Draw("HIST SAME");
    cdwp_tt_3pts_x_diff_per_layer[2]->Draw("HIST SAME");
    laplace_cdwp_tt_3pts_x_diff->SetLineWidth(3);
    laplace_cdwp_tt_3pts_x_diff->SetLineColor(kRed+1);
    laplace_cdwp_tt_3pts_x_diff->Draw("SAME");

    std::cout << "(CD/WP - TT) - x - 3 points resolution = " << laplace_cdwp_tt_3pts_x_diff->GetParameter(2) * std::sqrt(2.0) << " ± " << laplace_cdwp_tt_3pts_x_diff->GetParError(2) * std::sqrt(2.0) << std::endl;

    TLegend* leg_cdwp_tt_3pts_x_diff = new TLegend(0.15, 0.65, 0.35, 0.85);
    leg_cdwp_tt_3pts_x_diff->AddEntry(h_cdwp_tt_3pts_x_diff, "All", "l");
    leg_cdwp_tt_3pts_x_diff->AddEntry(cdwp_tt_3pts_x_diff_per_layer[0], "Layer 1", "l");
    leg_cdwp_tt_3pts_x_diff->AddEntry(cdwp_tt_3pts_x_diff_per_layer[1], "Layer 2", "l");
    leg_cdwp_tt_3pts_x_diff->AddEntry(cdwp_tt_3pts_x_diff_per_layer[2], "Layer 3", "l");
    leg_cdwp_tt_3pts_x_diff->AddEntry(laplace_cdwp_tt_3pts_x_diff, "Fit", "l");
    leg_cdwp_tt_3pts_x_diff->Draw();

    c_cdwp_tt_3pts_x_diff->SetTickx();
    c_cdwp_tt_3pts_x_diff->SetTicky();
    c_cdwp_tt_3pts_x_diff->SetLogy();
    c_cdwp_tt_3pts_x_diff->Update();

    TPaveStats* st_cdwp_tt_3pts_x_diff = (TPaveStats*)h_cdwp_tt_3pts_x_diff->FindObject("stats");
    st_cdwp_tt_3pts_x_diff->SetOptStat(ToROOTOpt(StatOpt::None));
    st_cdwp_tt_3pts_x_diff->SetOptFit(ToROOTOpt(FitOpt::Proba | FitOpt::Chi2NDF | FitOpt::AllParams | FitOpt::Errors));
    st_cdwp_tt_3pts_x_diff->SetX1NDC(0.55);
    st_cdwp_tt_3pts_x_diff->SetX2NDC(0.85);
    st_cdwp_tt_3pts_x_diff->SetY1NDC(0.65);
    st_cdwp_tt_3pts_x_diff->SetY2NDC(0.85);

    c_cdwp_tt_3pts_x_diff->Modified();
    c_cdwp_tt_3pts_x_diff->Update();

    // ============================================================================================
    // (CD/WP - TT) - y - 3 points
    // ============================================================================================

    TCanvas* c_cdwp_tt_3pts_y_diff = new TCanvas("c_cdwp_tt_3pts_y_diff", "c_cdwp_tt_3pts_y_diff", 1000, 1000);
    c_cdwp_tt_3pts_y_diff->cd();

    style(h_cdwp_tt_3pts_y_diff, kBlack, kSolid, 3);
    style(cdwp_tt_3pts_y_diff_per_layer[0], kBlue+1, kDashed, 3);
    style(cdwp_tt_3pts_y_diff_per_layer[1], kRed+1, kDashed, 3);
    style(cdwp_tt_3pts_y_diff_per_layer[2], kGreen-2, kDashed, 3);

    // TF1* laplace_cdwp_tt_3pts_y_diff = new TF1("laplace_cdwp_tt_3pts_y_diff","[0]*exp(-abs(x-[1])/[2]) + [3]", -5.0, 5.0);
    // laplace_cdwp_tt_3pts_y_diff->SetParameters(100.0, 0.0, 0.1, 1.0); // [0]=A, [1]=mu, [2]=tau, [3]=C
    TF1* laplace_cdwp_tt_3pts_y_diff = new TF1("laplace_cdwp_tt_3pts_y_diff","[0]*exp(-abs(x-[1])/[2]) + [3]*exp(-abs(x-[4])/[5])", -5.0, 5.0);
    laplace_cdwp_tt_3pts_y_diff->SetParameters(
        h_cdwp_tt_3pts_y_diff->GetMaximum(), 0.0, 0.5, 
        0.25 * h_cdwp_tt_3pts_y_diff->GetMaximum(), 0.0, 0.25
    );

    h_cdwp_tt_3pts_y_diff->Fit(laplace_cdwp_tt_3pts_y_diff, "R");
    h_cdwp_tt_3pts_y_diff->GetXaxis()->SetTitle("#Delta y (m) [CD/WP - TT]");
    h_cdwp_tt_3pts_y_diff->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_3pts_y_diff->GetYaxis()->SetTitle("Entries");
    h_cdwp_tt_3pts_y_diff->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_3pts_y_diff->GetYaxis()->SetTitleOffset(1.25);
    h_cdwp_tt_3pts_y_diff->Draw("HIST");
    cdwp_tt_3pts_y_diff_per_layer[0]->Draw("HIST SAME");
    cdwp_tt_3pts_y_diff_per_layer[1]->Draw("HIST SAME");
    cdwp_tt_3pts_y_diff_per_layer[2]->Draw("HIST SAME");
    laplace_cdwp_tt_3pts_y_diff->SetLineWidth(3);
    laplace_cdwp_tt_3pts_y_diff->SetLineColor(kRed+1);
    laplace_cdwp_tt_3pts_y_diff->Draw("SAME");
    
    std::cout << "(CD/WP - TT) - y - 3 points resolution = " << laplace_cdwp_tt_3pts_y_diff->GetParameter(2) * std::sqrt(2.0) << " ± " << laplace_cdwp_tt_3pts_y_diff->GetParError(2) * std::sqrt(2.0) << std::endl;

    TLegend* leg_cdwp_tt_3pts_y_diff = new TLegend(0.15, 0.65, 0.35, 0.85);
    leg_cdwp_tt_3pts_y_diff->AddEntry(h_cdwp_tt_3pts_y_diff, "All", "l");
    leg_cdwp_tt_3pts_y_diff->AddEntry(cdwp_tt_3pts_y_diff_per_layer[0], "Layer 1", "l");
    leg_cdwp_tt_3pts_y_diff->AddEntry(cdwp_tt_3pts_y_diff_per_layer[1], "Layer 2", "l");
    leg_cdwp_tt_3pts_y_diff->AddEntry(cdwp_tt_3pts_y_diff_per_layer[2], "Layer 3", "l");
    leg_cdwp_tt_3pts_y_diff->AddEntry(laplace_cdwp_tt_3pts_y_diff, "Fit", "l");
    leg_cdwp_tt_3pts_y_diff->Draw();

    c_cdwp_tt_3pts_y_diff->SetTickx();
    c_cdwp_tt_3pts_y_diff->SetTicky();
    c_cdwp_tt_3pts_y_diff->SetLogy();
    c_cdwp_tt_3pts_y_diff->Update();

     TPaveStats* st_cdwp_tt_3pts_y_diff = (TPaveStats*)h_cdwp_tt_3pts_y_diff->FindObject("stats");
    st_cdwp_tt_3pts_y_diff->SetOptStat(ToROOTOpt(StatOpt::None));
    st_cdwp_tt_3pts_y_diff->SetOptFit(ToROOTOpt(FitOpt::Proba | FitOpt::Chi2NDF | FitOpt::AllParams | FitOpt::Errors));
    st_cdwp_tt_3pts_y_diff->SetX1NDC(0.55);
    st_cdwp_tt_3pts_y_diff->SetX2NDC(0.85);
    st_cdwp_tt_3pts_y_diff->SetY1NDC(0.65);
    st_cdwp_tt_3pts_y_diff->SetY2NDC(0.85);

    c_cdwp_tt_3pts_y_diff->Modified();
    c_cdwp_tt_3pts_y_diff->Update();

    // ============================================================================================
    // (CD/WP - TT) - x - 2 points
    // ============================================================================================

    TCanvas* c_cdwp_tt_2pts_x_diff = new TCanvas("c_cdwp_tt_2pts_x_diff", "c_cdwp_tt_2pts_x_diff", 1000, 1000);
    c_cdwp_tt_2pts_x_diff->cd();

    style(h_cdwp_tt_2pts_x_diff, kBlack, kSolid, 3);
    style(cdwp_tt_2pts_x_diff_per_layer[0], kBlue+1, kDashed, 3);
    style(cdwp_tt_2pts_x_diff_per_layer[1], kRed+1, kDashed, 3);
    style(cdwp_tt_2pts_x_diff_per_layer[2], kGreen-2, kDashed, 3);

    // TF1* laplace_cdwp_tt_2pts_x_diff = new TF1("laplace_cdwp_tt_2pts_x_diff","[0]*exp(-abs(x-[1])/[2]) + [3]", -5.0, 5.0);
    // laplace_cdwp_tt_2pts_x_diff->SetParameters(100.0, 0.0, 0.1, 1.0); // [0]=A, [1]=mu, [2]=tau, [3]=C
    TF1* laplace_cdwp_tt_2pts_x_diff = new TF1("laplace_cdwp_tt_2pts_x_diff","[0]*exp(-abs(x-[1])/[2]) + [3]*exp(-abs(x-[4])/[5])", -5.0, 5.0);
    laplace_cdwp_tt_2pts_x_diff->SetParameters(
        h_cdwp_tt_2pts_x_diff->GetMaximum(), 0.0, 0.5, 
        0.25 * h_cdwp_tt_2pts_x_diff->GetMaximum(), 0.0, 0.25
    );

    h_cdwp_tt_2pts_x_diff->Fit(laplace_cdwp_tt_2pts_x_diff, "R");
    h_cdwp_tt_2pts_x_diff->GetXaxis()->SetTitle("#Delta x (m) [CD/WP - TT]");
    h_cdwp_tt_2pts_x_diff->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_2pts_x_diff->GetYaxis()->SetTitle("Entries");
    h_cdwp_tt_2pts_x_diff->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_2pts_x_diff->GetYaxis()->SetTitleOffset(1.25);
    h_cdwp_tt_2pts_x_diff->Draw("HIST");
    cdwp_tt_2pts_x_diff_per_layer[0]->Draw("HIST SAME");
    cdwp_tt_2pts_x_diff_per_layer[1]->Draw("HIST SAME");
    cdwp_tt_2pts_x_diff_per_layer[2]->Draw("HIST SAME");
    laplace_cdwp_tt_2pts_x_diff->SetLineWidth(3);
    laplace_cdwp_tt_2pts_x_diff->SetLineColor(kRed+1);
    laplace_cdwp_tt_2pts_x_diff->Draw("SAME");
    
    std::cout << "(CD/WP - TT) - x - 2 points resolution = " << laplace_cdwp_tt_2pts_x_diff->GetParameter(2) * std::sqrt(2.0) << " ± " << laplace_cdwp_tt_2pts_x_diff->GetParError(2) * std::sqrt(2.0) << std::endl;

    TLegend* leg_cdwp_tt_2pts_x_diff = new TLegend(0.15, 0.65, 0.35, 0.85);
    leg_cdwp_tt_2pts_x_diff->AddEntry(h_cdwp_tt_2pts_x_diff, "All", "l");
    leg_cdwp_tt_2pts_x_diff->AddEntry(cdwp_tt_2pts_x_diff_per_layer[0], "Layer 1", "l");
    leg_cdwp_tt_2pts_x_diff->AddEntry(cdwp_tt_2pts_x_diff_per_layer[1], "Layer 2", "l");
    leg_cdwp_tt_2pts_x_diff->AddEntry(cdwp_tt_2pts_x_diff_per_layer[2], "Layer 3", "l");
    leg_cdwp_tt_2pts_x_diff->AddEntry(laplace_cdwp_tt_2pts_x_diff, "Fit", "l");
    leg_cdwp_tt_2pts_x_diff->Draw();

    c_cdwp_tt_2pts_x_diff->SetTickx();
    c_cdwp_tt_2pts_x_diff->SetTicky();
    c_cdwp_tt_2pts_x_diff->SetLogy();
    c_cdwp_tt_2pts_x_diff->Update();

    TPaveStats* st_cdwp_tt_2pts_x_diff = (TPaveStats*)h_cdwp_tt_2pts_x_diff->FindObject("stats");
    st_cdwp_tt_2pts_x_diff->SetOptStat(ToROOTOpt(StatOpt::None));
    st_cdwp_tt_2pts_x_diff->SetOptFit(ToROOTOpt(FitOpt::Proba | FitOpt::Chi2NDF | FitOpt::AllParams | FitOpt::Errors));
    st_cdwp_tt_2pts_x_diff->SetX1NDC(0.55);
    st_cdwp_tt_2pts_x_diff->SetX2NDC(0.85);
    st_cdwp_tt_2pts_x_diff->SetY1NDC(0.65);
    st_cdwp_tt_2pts_x_diff->SetY2NDC(0.85);

    c_cdwp_tt_2pts_x_diff->Modified();
    c_cdwp_tt_2pts_x_diff->Update();

    // ============================================================================================
    // (CD/WP - TT) - y - 2 points
    // ============================================================================================

    TCanvas* c_cdwp_tt_2pts_y_diff = new TCanvas("c_cdwp_tt_2pts_y_diff", "c_cdwp_tt_2pts_y_diff", 1000, 1000);
    c_cdwp_tt_2pts_y_diff->cd();

    style(h_cdwp_tt_2pts_y_diff, kBlack, kSolid, 3);
    style(cdwp_tt_2pts_y_diff_per_layer[0], kBlue+1, kDashed, 3);
    style(cdwp_tt_2pts_y_diff_per_layer[1], kRed+1, kDashed, 3);
    style(cdwp_tt_2pts_y_diff_per_layer[2], kGreen-2, kDashed, 3);

    // TF1* laplace_cdwp_tt_2pts_y_diff = new TF1("laplace_cdwp_tt_2pts_y_diff","[0]*exp(-abs(x-[1])/[2]) + [3]", -5.0, 5.0);
    // laplace_cdwp_tt_2pts_y_diff->SetParameters(100.0, 0.0, 0.1, 1.0); // [0]=A, [1]=mu, [2]=tau, [3]=C
    TF1* laplace_cdwp_tt_2pts_y_diff = new TF1("laplace_cdwp_tt_2pts_y_diff","[0]*exp(-abs(x-[1])/[2]) + [3]*exp(-abs(x-[4])/[5])", -5.0, 5.0);
    laplace_cdwp_tt_2pts_y_diff->SetParameters(
        h_cdwp_tt_2pts_y_diff->GetMaximum(), 0.0, 0.5, 
        0.25 * h_cdwp_tt_2pts_y_diff->GetMaximum(), 0.0, 0.25
    );

    h_cdwp_tt_2pts_y_diff->Fit(laplace_cdwp_tt_2pts_y_diff, "R");
    h_cdwp_tt_2pts_y_diff->GetXaxis()->SetTitle("#Delta y (m) [CD/WP - TT]");
    h_cdwp_tt_2pts_y_diff->GetXaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_2pts_y_diff->GetYaxis()->SetTitle("Entries");
    h_cdwp_tt_2pts_y_diff->GetYaxis()->CenterTitle(kTRUE);
    h_cdwp_tt_2pts_y_diff->GetYaxis()->SetTitleOffset(1.25);
    h_cdwp_tt_2pts_y_diff->Draw("HIST");
    cdwp_tt_2pts_y_diff_per_layer[0]->Draw("HIST SAME");
    cdwp_tt_2pts_y_diff_per_layer[1]->Draw("HIST SAME");
    cdwp_tt_2pts_y_diff_per_layer[2]->Draw("HIST SAME");
    laplace_cdwp_tt_2pts_y_diff->SetLineWidth(3);
    laplace_cdwp_tt_2pts_y_diff->SetLineColor(kRed+1);
    laplace_cdwp_tt_2pts_y_diff->Draw("SAME");

    std::cout << "(CD/WP - TT) - y - 2 points resolution = " << laplace_cdwp_tt_2pts_y_diff->GetParameter(2) * std::sqrt(2.0) << " ± " << laplace_cdwp_tt_2pts_y_diff->GetParError(2) * std::sqrt(2.0) << std::endl;

    TLegend* leg_cdwp_tt_2pts_y_diff = new TLegend(0.15, 0.65, 0.35, 0.85);
    leg_cdwp_tt_2pts_y_diff->AddEntry(h_cdwp_tt_2pts_y_diff, "All", "l");
    leg_cdwp_tt_2pts_y_diff->AddEntry(cdwp_tt_2pts_y_diff_per_layer[0], "Layer 1", "l");
    leg_cdwp_tt_2pts_y_diff->AddEntry(cdwp_tt_2pts_y_diff_per_layer[1], "Layer 2", "l");
    leg_cdwp_tt_2pts_y_diff->AddEntry(cdwp_tt_2pts_y_diff_per_layer[2], "Layer 3", "l");
    leg_cdwp_tt_2pts_y_diff->AddEntry(laplace_cdwp_tt_2pts_y_diff, "Fit", "l");
    leg_cdwp_tt_2pts_y_diff->Draw();

    c_cdwp_tt_2pts_y_diff->SetTickx();
    c_cdwp_tt_2pts_y_diff->SetTicky();
    c_cdwp_tt_2pts_y_diff->SetLogy();
    c_cdwp_tt_2pts_y_diff->Update();

    TPaveStats* st_cdwp_tt_2pts_y_diff = (TPaveStats*)h_cdwp_tt_2pts_y_diff->FindObject("stats");
    st_cdwp_tt_2pts_y_diff->SetOptStat(ToROOTOpt(StatOpt::None));
    st_cdwp_tt_2pts_y_diff->SetOptFit(ToROOTOpt(FitOpt::Proba | FitOpt::Chi2NDF | FitOpt::AllParams | FitOpt::Errors));
    st_cdwp_tt_2pts_y_diff->SetX1NDC(0.55);
    st_cdwp_tt_2pts_y_diff->SetX2NDC(0.85);
    st_cdwp_tt_2pts_y_diff->SetY1NDC(0.65);
    st_cdwp_tt_2pts_y_diff->SetY2NDC(0.85);

    c_cdwp_tt_2pts_y_diff->Modified();
    c_cdwp_tt_2pts_y_diff->Update();

    // ============================================================================================
    // Missing layers
    // ============================================================================================

    TCanvas* c_missing_layer = new TCanvas("c_missing_layer", "c_missing_layer", 1000, 1000);
    c_missing_layer->cd();

    style(h_missing_layer, kBlue+1, kSolid, 3);

    h_missing_layer->SetStats(0);
    h_missing_layer->GetXaxis()->SetTitle("Missing layer");
    h_missing_layer->GetXaxis()->CenterTitle(kTRUE);
    h_missing_layer->GetYaxis()->SetTitle("Entries");
    h_missing_layer->GetYaxis()->CenterTitle(kTRUE);
    h_missing_layer->GetYaxis()->SetTitleOffset(1.5);
    h_missing_layer->Draw();

    std::cout << "Missing layers: | ";
    for (int i = 1; i <= 6; ++i) {
        std::cout << i << " = " << h_missing_layer->GetBinContent(i) << " | ";
    }
    std::cout << '\n';

    c_missing_layer->SetTickx();
    c_missing_layer->SetTicky();
    c_missing_layer->Update();

    // ============================================================================================
    // X-Y map for each main layer - No TT
    // ============================================================================================

    TCanvas* c_xy_cdwp_no_tt_stack = new TCanvas("c_xy_cdwp_no_tt_stack", "XY Points per Layer - No TT", 800, 1000);
    c_xy_cdwp_no_tt_stack->Divide(1, 3);

    int layers_no_tt[3] = {2, 1, 0};

    for (int p = 0; p < 3; p++) {
        TVirtualPad* pad = c_xy_cdwp_no_tt_stack->cd(p + 1);
        pad->SetBottomMargin(0.20);

        TH2D* h = xy_pts_per_layer_cdwp_no_tt[layers_no_tt[p]];
        h->SetStats(0);
        h->SetTitle("");
        h->GetXaxis()->SetLabelSize(0.07);
        h->GetYaxis()->SetLabelSize(0.07);
        h->GetZaxis()->SetLabelSize(0.07);
        h->GetXaxis()->SetTitleSize(0.07);
        h->GetYaxis()->SetTitleSize(0.07);
        h->GetZaxis()->SetTitleSize(0.07);
        h->GetXaxis()->SetTitleOffset(1.0);
        h->GetYaxis()->SetTitleOffset(0.5);
        h->GetZaxis()->SetTitleOffset(0.25);
        h->GetXaxis()->CenterTitle(kTRUE);
        h->GetYaxis()->CenterTitle(kTRUE);
        h->GetXaxis()->SetTitle("x (m)");
        h->GetYaxis()->SetTitle("y (m)");
        h->GetZaxis()->SetTitle("Entries");
        h->Draw("COLZ");

        TBox* rect_outer = new TBox(-23.5, -10.0, 23.5, 10.0);
        rect_outer->SetFillStyle(0);
        rect_outer->SetLineColor(kBlack);
        rect_outer->SetLineWidth(2);
        rect_outer->SetLineStyle(kDotted);
        rect_outer->Draw("SAME");

        TBox* rect_inner = new TBox(-3.35, -3.35, 3.35, 3.35);
        rect_inner->SetFillStyle(0);
        rect_inner->SetLineColor(kBlack);
        rect_inner->SetLineWidth(2);
        rect_inner->SetLineStyle(kDotted);
        rect_inner->Draw("SAME");

        pad->SetTickx();
        pad->SetTicky();
    }

    // ============================================================================================
    // X-Y map for each main layer - No TT
    // ============================================================================================

    TCanvas* c_xy_cdwp_no_tt_zoom_stack = new TCanvas("c_xy_cdwp_no_tt_zoom_stack", "Zoomed XY Points per Layer - No TT", 800, 1000);
    c_xy_cdwp_no_tt_zoom_stack->Divide(1, 3);

    for (int p = 0; p < 3; p++) {
        TVirtualPad* pad = c_xy_cdwp_no_tt_zoom_stack->cd(p + 1);
        pad->SetBottomMargin(0.20);

        TH2D* h = xy_pts_per_layer_cdwp_no_tt_zoom[layers_no_tt[p]];
        h->SetStats(0);
        h->SetTitle("");
        h->GetXaxis()->SetLabelSize(0.07);
        h->GetYaxis()->SetLabelSize(0.07);
        h->GetZaxis()->SetLabelSize(0.07);
        h->GetXaxis()->SetTitleSize(0.07);
        h->GetYaxis()->SetTitleSize(0.07);
        h->GetZaxis()->SetTitleSize(0.07);
        h->GetXaxis()->SetTitleOffset(1.0);
        h->GetYaxis()->SetTitleOffset(0.5);
        h->GetZaxis()->SetTitleOffset(0.25);
        h->GetXaxis()->CenterTitle(kTRUE);
        h->GetYaxis()->CenterTitle(kTRUE);
        h->GetXaxis()->SetTitle("x (m)");
        h->GetYaxis()->SetTitle("y (m)");
        h->GetZaxis()->SetTitle("Entries");
        h->Draw("COLZ");

        pad->SetTickx();
        pad->SetTicky();
    }

    // ============================================================================================
    // X-Y map for each main layer - No TT
    // ============================================================================================

    TCanvas* c_xy_cdwp_2pts_extrapolation = new TCanvas("c_xy_cdwp_2pts_extrapolation", "Zoomed XY Points per Layer - 2 points TT", 800, 1000);
    c_xy_cdwp_2pts_extrapolation->Divide(1, 3);

    for (int p = 0; p < 3; p++) {
        TVirtualPad* pad = c_xy_cdwp_2pts_extrapolation->cd(p + 1);
        pad->SetBottomMargin(0.20);

        TH2D* h = xy_pts_per_layer_cdwp_2pts_extrapolation[layers_no_tt[p]];

        h->SetStats(0);
        h->SetTitle("");
        h->GetXaxis()->SetLabelSize(0.07);
        h->GetYaxis()->SetLabelSize(0.07);
        h->GetZaxis()->SetLabelSize(0.07);
        h->GetXaxis()->SetTitleSize(0.07);
        h->GetYaxis()->SetTitleSize(0.07);
        h->GetZaxis()->SetTitleSize(0.07);
        h->GetXaxis()->SetTitleOffset(1.0);
        h->GetYaxis()->SetTitleOffset(0.5);
        h->GetZaxis()->SetTitleOffset(0.25);
        h->GetXaxis()->CenterTitle(kTRUE);
        h->GetYaxis()->CenterTitle(kTRUE);
        h->GetXaxis()->SetTitle("x (m)");
        h->GetYaxis()->SetTitle("y (m)");
        h->GetZaxis()->SetTitle("Entries");
        h->Draw("COLZ");

        pad->SetTickx();
        pad->SetTicky();
    }
}