#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <unordered_set>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TVector3.h>

inline TTimeStamp operator-(const TTimeStamp& lhs, const TTimeStamp& rhs) {
    return TTimeStamp(lhs.GetSec() - rhs.GetSec(), lhs.GetNanoSec() - rhs.GetNanoSec());
}

inline TTimeStamp operator+(const TTimeStamp& lhs, const TTimeStamp& rhs) {
    return TTimeStamp(lhs.GetSec() + rhs.GetSec(), lhs.GetNanoSec() + rhs.GetNanoSec());
}

void extract_plot_from_reco_matches(const char* filename) {
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

    double cdwp_tt_coordinate_offset = 26452.0; // 26470 or 26452
    double layer_0_zmin = 24000.0, layer_0_zmax = 25000.0; // [24000, 25000] 
    double layer_1_zmin = 25500.0, layer_1_zmax = 26500.0; // [25500, 26500] 
    double layer_2_zmin = 27000.0, layer_2_zmax = 28000.0; // [27000, 28000] 
    double layer_3_zmin = 30000.0, layer_3_zmax = 30200.0; // [30000, 30200] 
    double layer_4_zmin = 30200.0, layer_4_zmax = 30400.0; // [30200, 30400] 
    double layer_5_zmin = 30400.0, layer_5_zmax = 30600.0; // [30400, 30600] 

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

    std::vector<int> run_ids;
    std::vector<double> angles;
    std::vector<double> distances;
    std::vector<double> clippingness;
    std::vector<double> chi2_cdwp_v;
    std::vector<unsigned char> det_cdwp_v;
    std::vector<double> chi2_tt_v;
    std::vector<double> zenith_cdwp_v;
    std::vector<double> azimuth_cdwp_v;
    std::vector<double> zenith_tt_v;
    std::vector<double> azimuth_tt_v;

    TVector3 ipos, fpos;
    TVector3 pos_cdwp, dir_cdwp, pos_tt, dir_tt;

    std::cout << "[INFO] Number of entries: " << tree->GetEntries() << '\n';

    TH1D* h_ts_diff = new TH1D("h_ts_diff", "h_ts_diff", 100, -1000.0, 1000.0);
    TH1D* h_z_pt = new TH1D("h_z_pt", "h_z_pt", 1000, -10000.0, 10000.0 /* 24000.0, 31000.0 */);

    auto layer_id = [&](double z) {
        if (layer_0_zmin <= z && z <= layer_0_zmax) return 0;  // main
        if (layer_1_zmin <= z && z <= layer_1_zmax) return 1;  // main
        if (layer_2_zmin <= z && z <= layer_2_zmax) return 2;  // main
        if (layer_3_zmin <= z && z <= layer_3_zmax) return 3;  // chimney
        if (layer_4_zmin <= z && z <= layer_4_zmax) return 4;  // chimney
        if (layer_5_zmin <= z && z <= layer_5_zmax) return 5;  // chimney

        return -1; // not inside any valid layer
    };

    for (int k = 0; k < tree->GetEntries(); ++k) {
        if (k % 1000 == 0) std::cout << "\rEntries: " << k << " / " << tree->GetEntries();
        tree->GetEntry(k);

        if (NTotPoints == 0) continue;
        if (NTracks != 1) continue;
        if (NPoints[0] < 3) continue;

        std::unordered_set<int> layers_hit;
        layers_hit.reserve(6);

        h_ts_diff->Fill(1.0e9 * (TTimeStamp{sec, nsec} - *start_TS));
        for (int i = 0; i < NTotPoints; ++i) {
            h_z_pt->Fill(PointZ[i] + cdwp_tt_coordinate_offset);
            int lid = layer_id(PointZ[i] + cdwp_tt_coordinate_offset);
            if (lid >= 0) layers_hit.insert(lid);
        }

        if (layers_hit.size() < 3) continue; // at least 3 different layers

        pos_tt.SetXYZ(Coeff0[0], Coeff1[0], Coeff2[0] + cdwp_tt_coordinate_offset);
        dir_tt.SetXYZ(Coeff3[0], Coeff4[0], Coeff5[0]);
        dir_tt = dir_tt.Unit();

        if (dir_tt.Cross(pos_tt).Mag() / 1000.0 > 17.7) continue;

        std::unordered_map<std::string, std::size_t> method_count = {
            {"CdWpTtChi2", 0ul},
            {"WpClassify", 0ul}
        };
        std::size_t j = 0ul;
        for (std::size_t i = 0ul; i < method->size(); ++i) {
            ++(method_count[(*method)[i]]);
            if ((*method)[i] == "CdWpTtChi2") j = i;
        }
        if (method_count["WpClassify"] != 1ul || method_count["CdWpTtChi2"] != 1ul) continue;

        ipos.SetXYZ((*iposx)[j], (*iposy)[j], (*iposz)[j]);
        fpos.SetXYZ((*fposx)[j], (*fposy)[j], (*fposz)[j]);
        pos_cdwp = ipos;
        dir_cdwp = (fpos - ipos).Unit();

        run_ids.push_back(run_id);
        angles.push_back(dir_tt.Angle(dir_cdwp) * 180.0 / M_PI);
        distances.push_back(
            ((pos_tt - (pos_tt * dir_tt) * dir_tt) - (pos_cdwp - (pos_cdwp * dir_cdwp) * dir_cdwp)).Mag() / 1000.0
        );
        clippingness.push_back(
            dir_tt.Cross(pos_tt).Mag() / 1000.0
        );
        chi2_cdwp_v.push_back((*quality)[j]);
        chi2_tt_v.push_back(Chi2[0]);
        det_cdwp_v.push_back((*det)[j]);

        zenith_cdwp_v.push_back(-dir_cdwp.CosTheta());
        azimuth_cdwp_v.push_back( 180.0 / M_PI * (dir_cdwp.Phi() > 0.0 ? dir_cdwp.Phi() : dir_cdwp.Phi() + 2.0 * M_PI) );
        zenith_tt_v.push_back(-dir_tt.CosTheta());
        azimuth_tt_v.push_back( 180.0 / M_PI * (dir_tt.Phi() > 0.0 ? dir_tt.Phi() : dir_tt.Phi() + 2.0 * M_PI) );

        // if (angles.back() > 3.0 || distances.back() > 1.0) {
        //     std::cout << TTimeStamp{sec, nsec} << ", angle: " << angles.back() << ", distance: " << distances.back() 
        //               << ", cdwp_dir: (" << dir_cdwp.Theta() << ", " << dir_cdwp.Phi() << "), tt_dir: (" << dir_tt.Theta() << ", "<< dir_tt.Phi() << ")\n";
            // std::cout << cdwp_sec_out << ' ' << cdwp_nsec_out << ' ' 
            //           << tt_x << ' ' << tt_y << ' ' << tt_z << ' ' << tt_dx << ' ' << tt_dy << ' ' << tt_dz << ' '
            //           << cdwp_x << ' ' << cdwp_y << ' ' << cdwp_z << ' ' << cdwp_dx << ' ' << cdwp_dy << ' ' << cdwp_dz << '\n';
            // std::cout << cdwp_sec_out << ' ' << cdwp_nsec_out << ' ' << tt_x << ' ' << tt_y << ' ' << tt_z << ' ' << tt_dx << ' ' << tt_dy << ' ' << tt_dz << '\n';
        // }
    }

    TCanvas* c_ts_diff = new TCanvas("c_ts_diff", "c_ts_diff", 1000, 1000);
    c_ts_diff->cd();
    h_ts_diff->Draw();
    c_ts_diff->Update();

    TCanvas* c_z_pt = new TCanvas("c_z_pt", "c_z_pt", 1000, 1000);
    c_z_pt->cd();
    h_z_pt->Draw();
    c_z_pt->Update();

    TH1I* h_det = new TH1I("h_det", "h_det", 8, 0, 8);
    TH1D* h_angle = new TH1D("h_angle", "h_angle", 20, 0.0, 5.0);
    TH1D* h_distance = new TH1D("h_distance", "h_distance", 20, 0.0, 2.0);
    TH2D* h_angle_clippingness = new TH2D("h_angle_clippingness", "h_angle_clippingness", 50, 0.0, 90.0, 50, 0.0, 20.5); // [0.0, 20.0] or [0.0, 90.0] or [0.0, 180.0]
    TH2D* h_distance_clippingness = new TH2D("h_distance_clippingness", "h_distance_clippingness", 50, 0.0, 20.0, 50, 0.0, 20.5); // [0.0, 5.0] or [0.0, 20.0] or [0.0, 40.0]
    TH2D* h_angle_distance = new TH2D("h_angle_distance", "h_angle_distance", 50, 0.0, 90.0, 50, 0.0, 20.0);
    TH2D* h_angle_chi2cdwp = new TH2D("h_angle_chi2cdwp", "h_angle_chi2cdwp", 50, 0.0, 20.0, 50, 0.0, 10.0);
    TH2D* h_distance_chi2cdwp = new TH2D("h_distance_chi2cdwp", "h_distance_chi2cdwp", 50, 0.0, 5.0, 50, 0.0, 10.0);

    // std::unordered_map<int, TH1D*> map_angle_run;
    // std::unordered_map<int, TH1D*> map_distance_run;
    TH1D* h_angle_run_summry = new TH1D("h_angle_run_summary", "h_angle_run_summary", max_run - min_run + 1, min_run - 0.5, max_run + 0.5);
    TH1D* h_distance_run_summry = new TH1D("h_distance_run_summary", "h_distance_run_summary", max_run - min_run + 1, min_run - 0.5, max_run + 0.5);
    std::unordered_map<int, std::vector<double>> map_angle_run_values;
    std::unordered_map<int, std::vector<double>> map_distance_run_values;
    std::unordered_map<int, int> map_angle_run_counts;
    std::unordered_map<int, int> map_distance_run_counts;

    for (std::size_t k = 0ul; k < angles.size(); ++k) {
        // if (chi2_cdwp_v[k] > 0.0) continue;
        h_det->Fill(det_cdwp_v[k]);
        h_angle->Fill(angles[k]);
        h_distance->Fill(distances[k]);
        h_angle_clippingness->Fill(angles[k], clippingness[k]);
        h_distance_clippingness->Fill(distances[k], clippingness[k]);
        h_angle_distance->Fill(angles[k], distances[k]);
        h_angle_chi2cdwp->Fill(angles[k], chi2_cdwp_v[k]);
        h_distance_chi2cdwp->Fill(distances[k], chi2_cdwp_v[k]);

        // if (map_angle_run.find(run_ids[k]) == map_angle_run.end()) {
        //     map_angle_run[run_ids[k]] = new TH1D(Form("h_angle_run_%d", run_ids[k]), Form("h_angle_run_%d", run_ids[k]), 20, 0.0, 5.0);
        // }
        // map_angle_run[run_ids[k]]->Fill(angles[k]);
        // if (map_distance_run.find(run_ids[k]) == map_distance_run.end()) {
        //     map_distance_run[run_ids[k]] = new TH1D(Form("h_distance_run_%d", run_ids[k]), Form("h_distance_run_%d", run_ids[k]), 20, 0.0, 2.0);
        // }
        // map_distance_run[run_ids[k]]->Fill(distances[k]);
        map_angle_run_values[run_ids[k]].push_back(angles[k]);
        map_distance_run_values[run_ids[k]].push_back(distances[k]);
        ++(map_angle_run_counts[run_ids[k]]);
        ++(map_distance_run_counts[run_ids[k]]);
    }
    for (const auto& [run_id, values] : map_angle_run_values) {
        std::nth_element(values.begin(), values.begin() + values.size() * 682 / 1000, values.end());
        h_angle_run_summry->SetBinContent(run_id - min_run + 1, angles[angles.size() * 682 / 1000]);
    }
    for (const auto& [run_id, values] : map_distance_run_values) {
        std::nth_element(values.begin(), values.begin() + values.size() * 682 / 1000, values.end());
        h_distance_run_summry->SetBinContent(run_id - min_run + 1, distances[distances.size() * 682 / 1000]);
    }

    std::cout << h_angle->Integral(0, 20) << '/' << h_angle->Integral(0, 50) << '\n';

    h_angle->Scale(1.0 / h_angle->Integral());
    h_distance->Scale(1.0 / h_distance->Integral());

    double bottom_pixel_angle = 773.0;
    double top_pixel_angle = 175.0;
    double axis_scale_angle = 0.8;
    std::vector<double> sftm_angle_pixel =      {598.0, 362.0, 284.0, 352.0, 439.0, 529.0, 602.0, 625.0, 685.0, 700.0, 710.0, 719.0, 725.0, 725.0, 735.0, 736.0, 732.0, 730.0, 737.0, 735.0};
    std::vector<double> ml_angle_pixel =        {566.0, 235.0, 161.0, 322.0, 439.0, 579.0, 654.0, 697.0, 726.0, 748.0, 751.0, 735.0, 755.0, 737.0, 761.0, 751.0, 751.0, 758.0, 762.0, 752.0};
    std::vector<double> wpcluster_angle_pixel = {718.0, 585.0, 510.0, 484.0, 492.0, 501.0, 542.0, 589.0, 644.0, 624.0, 677.0, 688.0, 686.0, 696.0, 722.0, 722.0, 720.0, 741.0, 749.0, 735.0};
    std::vector<double> cdcluster_angle_pixel = {696.0, 535.0, 443.0, 381.0, 421.0, 466.0, 488.0, 571.0, 643.0, 665.0, 696.0, 723.0, 733.0, 743.0, 732.0, 744.0, 739.0, 750.0, 757.0, 758.0};
    double bottom_pixel_distance = 766.0;
    double top_pixel_distance = 143.0;
    double axis_scale_distance = 2.5;
    std::vector<double> sftm_distance_pixel =      {702.0, 484.0, 358.0, 324.0, 456.0, 567.0, 623.0, 672.0, 681.0, 699.0, 687.0, 702.0, 718.0, 730.0, 723.0, 735.0, 727.0, 734.0, 742.0, 746.0};
    std::vector<double> ml_distance_pixel =        {589.0, 256.0, 184.0, 375.0, 527.0, 602.0, 678.0, 700.0, 709.0, 726.0, 741.0, 755.0, 748.0, 735.0, 752.0, 747.0, 747.0, 752.0, 743.0, 743.0};
    std::vector<double> wpcluster_distance_pixel = {745.0, 618.0, 547.0, 517.0, 456.0, 503.0, 564.0, 609.0, 626.0, 611.0, 644.0, 650.0, 684.0, 687.0, 714.0, 726.0, 717.0, 732.0, 730.0, 734.0};
    std::vector<double> cdcluster_distance_pixel = {714.0, 469.0, 320.0, 333.0, 354.0, 545.0, 573.0, 681.0, 679.0, 730.0, 718.0, 742.0, 731.0, 738.0, 751.0, 746.0, 748.0, 751.0, 744.0, 747.0};

    TH1D* h_sftm_angle = new TH1D("h_sftm_angle", "h_sftm_angle", 20, 0.0, 5.0);
    TH1D* h_ml_angle = new TH1D("h_ml_angle", "h_ml_angle", 20, 0.0, 5.0); 
    TH1D* h_wpcluster_angle = new TH1D("h_wpcluster_angle", "h_wpcluster_angle", 20, 0.0, 5.0);
    TH1D* h_cdcluster_angle = new TH1D("h_cdcluster_angle", "h_cdcluster_angle", 20, 0.0, 5.0);
    TH1D* h_sftm_distance = new TH1D("h_sftm_distance", "h_sftm_distance", 20, 0.0, 2.0);
    TH1D* h_ml_distance = new TH1D("h_ml_distance", "h_ml_distance", 20, 0.0, 2.0); 
    TH1D* h_wpcluster_distance = new TH1D("h_wpcluster_distance", "h_wpcluster_distance", 20, 0.0, 2.0);
    TH1D* h_cdcluster_distance = new TH1D("h_cdcluster_distance", "h_cdcluster_distance", 20, 0.0, 2.0);
    for (int k = 0; k < 20; ++k) {
        h_sftm_angle->SetBinContent(k + 1, axis_scale_angle * (bottom_pixel_angle - sftm_angle_pixel[k]) / (bottom_pixel_angle - top_pixel_angle));
        h_ml_angle->SetBinContent(k + 1, axis_scale_angle * (bottom_pixel_angle - ml_angle_pixel[k]) / (bottom_pixel_angle - top_pixel_angle));
        h_wpcluster_angle->SetBinContent(k + 1, axis_scale_angle * (bottom_pixel_angle - wpcluster_angle_pixel[k]) / (bottom_pixel_angle - top_pixel_angle));
        h_cdcluster_angle->SetBinContent(k + 1, axis_scale_angle * (bottom_pixel_angle - cdcluster_angle_pixel[k]) / (bottom_pixel_angle - top_pixel_angle));
        h_sftm_distance->SetBinContent(k + 1, axis_scale_distance * (bottom_pixel_distance - sftm_distance_pixel[k]) / (bottom_pixel_distance - top_pixel_distance));
        h_ml_distance->SetBinContent(k + 1, axis_scale_distance * (bottom_pixel_distance - ml_distance_pixel[k]) / (bottom_pixel_distance - top_pixel_distance));
        h_wpcluster_distance->SetBinContent(k + 1, axis_scale_distance * (bottom_pixel_distance - wpcluster_distance_pixel[k]) / (bottom_pixel_distance - top_pixel_distance));
        h_cdcluster_distance->SetBinContent(k + 1, axis_scale_distance * (bottom_pixel_distance - cdcluster_distance_pixel[k]) / (bottom_pixel_distance - top_pixel_distance));
    }
    h_sftm_angle->Scale(1.0 / h_sftm_angle->Integral());
    h_ml_angle->Scale(1.0 / h_ml_angle->Integral());
    h_wpcluster_angle->Scale(1.0 / h_wpcluster_angle->Integral());
    h_cdcluster_angle->Scale(1.0 / h_cdcluster_angle->Integral());
    h_sftm_distance->Scale(1.0 / h_sftm_distance->Integral());
    h_ml_distance->Scale(1.0 / h_ml_distance->Integral());
    h_wpcluster_distance->Scale(1.0 / h_wpcluster_distance->Integral());
    h_cdcluster_distance->Scale(1.0 / h_cdcluster_distance->Integral());

    /* std::function<void(TH1*, Color_t, Style_t, int)> */ auto set_style = [](TH1* h, Color_t color, Style_t style, int alpha = 100) {
        h->SetLineColor(color);
        // h->SetFillColor(color + (alpha << 24));
        h->SetFillStyle(3004);       // transparent fill pattern
        h->SetLineWidth(4);
        h->SetLineStyle(style);
    };

    int min_run = *std::min_element(run_ids.begin(), run_ids.end());
    int max_run = *std::max_element(run_ids.begin(), run_ids.end());

    TCanvas* c_angle_run = new TCanvas("c_angle_run", "c_angle_run", 1000, 1000);
    c_angle_run->cd();
    // TH1D* h_angle_run_summry = new TH1D("h_angle_run_summary", "h_angle_run_summary", max_run - min_run + 1, min_run - 0.5, max_run + 0.5);
    // for (const auto& [run, h_run] : map_angle_run) {
    //     double probs[1] = {0.68};
    //     double quantiles[1];
    //     h_run->GetQuantiles(1, quantiles, probs);
    //     h_angle_run_summry->SetBinContent(run - min_run + 1, quantiles[0]);
    // }
    h_angle_run_summry->GetXaxis()->SetTitle("Run ID");
    h_angle_run_summry->GetXaxis()->CenterTitle(kTRUE);
    h_angle_run_summry->GetYaxis()->SetTitle("#alpha (deg) at 68% quantile");
    h_angle_run_summry->GetYaxis()->CenterTitle(kTRUE);
    h_angle_run_summry->GetYaxis()->SetTitleOffset(1.5);
    h_angle_run_summry->SetStats(0);
    h_angle_run_summry->Draw("HIST");
    c_angle_run->SetTickx();
    c_angle_run->SetTicky();
    c_angle_run->Update();

    TCanvas* c_distance_run = new TCanvas("c_distance_run", "c_distance_run", 1000, 1000);
    c_distance_run->cd();
    // TH1D* h_distance_run_summry = new TH1D("h_distance_run_summary", "h_distance_run_summary", max_run - min_run + 1, min_run - 0.5, max_run + 0.5);
    // for (const auto& [run, h_run] : map_distance_run) {
    //     double probs[1] = {0.68};
    //     double quantiles[1];
    //     h_run->GetQuantiles(1, quantiles, probs);
    //     h_distance_run_summry->SetBinContent(run - min_run + 1, quantiles[0]);
    // }
    h_distance_run_summry->GetXaxis()->SetTitle("Run ID");
    h_distance_run_summry->GetXaxis()->CenterTitle(kTRUE);
    h_distance_run_summry->GetYaxis()->SetTitle("d_{mid} (m) at 68% quantile");
    h_distance_run_summry->GetYaxis()->CenterTitle(kTRUE);
    h_distance_run_summry->GetYaxis()->SetTitleOffset(1.5);
    h_distance_run_summry->SetStats(0);
    h_distance_run_summry->Draw("HIST");
    c_distance_run->SetTickx();
    c_distance_run->SetTicky();
    c_distance_run->Update();

    TCanvas* c_det = new TCanvas("c_det", "c_det", 1000, 1000);
    c_det->cd();
    h_det->Draw();
    c_det->Update();

    TCanvas* c_angle = new TCanvas("c_angle", "c_angle", 1000, 1000);
    c_angle->cd();
    
    h_angle->SetStats(0);
    h_angle->SetMaximum(1.1 * std::max({h_angle->GetMaximum(), h_sftm_angle->GetMaximum(), h_ml_angle->GetMaximum()}));
    
    set_style(h_angle, kBlack, kSolid);
    set_style(h_sftm_angle, kOrange+2, kSolid);
    set_style(h_ml_angle, kRed, kDotted);
    set_style(h_wpcluster_angle, kViolet, kDashed);
    set_style(h_cdcluster_angle, kGreen+2, kDashed);

    h_angle->GetXaxis()->SetTitle("#alpha (deg)");
    h_angle->GetXaxis()->CenterTitle(kTRUE);
    h_angle->GetYaxis()->SetTitle("Entries");
    h_angle->GetYaxis()->CenterTitle(kTRUE);
    h_angle->GetYaxis()->SetTitleOffset(1.5);
    h_angle->Draw("HIST");
    h_sftm_angle->Draw("HIST SAME");
    h_ml_angle->Draw("HIST SAME");
    h_wpcluster_angle->Draw("HIST SAME");
    h_cdcluster_angle->Draw("HIST SAME");

    TLegend* leg_angle = new TLegend(0.55, 0.65, 0.85, 0.85);
    leg_angle->AddEntry(h_angle, "Joint #chi^2", "l");
    leg_angle->AddEntry(h_sftm_angle, "SFTM", "l");
    leg_angle->AddEntry(h_ml_angle, "ML", "l");
    leg_angle->AddEntry(h_wpcluster_angle, "WP cluster", "l");
    leg_angle->AddEntry(h_cdcluster_angle, "CD cluster", "l");
    leg_angle->Draw();

    double x_text = 0.55;
    double y_text = 0.63;   // Slightly below the legend's lower Y (0.65)

    TLatex* tex1_angle = new TLatex();
    tex1_angle->SetNDC();
    tex1_angle->SetTextFont(72);
    tex1_angle->SetTextSize(0.035);
    tex1_angle->SetTextAlign(13); // left-bottom alignment
    tex1_angle->DrawLatex(x_text, y_text, Form("RUN %d-%d", min_run, max_run));

    TLatex* tex2_angle = new TLatex();
    tex2_angle->SetNDC();
    tex2_angle->SetTextFont(42);
    tex2_angle->SetTextSize(0.025);
    tex2_angle->SetTextAlign(13); // left-bottom alignment
    tex2_angle->DrawLatex(x_text, y_text - 0.04, Form("Entries %lld", static_cast<long long>(h_angle->GetEntries())));

    c_angle->SetTickx();
    c_angle->SetTicky();
    c_angle->Update();

    TCanvas* c_distance = new TCanvas("c_distance", "c_distance", 1000, 1000);
    c_distance->cd();
    
    h_distance->SetStats(0);
    h_distance->SetMaximum(1.1 * std::max({h_angle->GetMaximum(), h_sftm_angle->GetMaximum(), h_ml_angle->GetMaximum()}));
    
    set_style(h_distance, kBlack, kSolid);
    set_style(h_sftm_distance, kOrange+2, kSolid);
    set_style(h_ml_distance, kRed, kDotted);
    set_style(h_wpcluster_distance, kViolet, kDashed);
    set_style(h_cdcluster_distance, kGreen+2, kDashed);
    
    h_distance->GetXaxis()->SetTitle("d_{mid} (m)");
    h_distance->GetXaxis()->CenterTitle(kTRUE);
    h_distance->GetYaxis()->SetTitle("Entries");
    h_distance->GetYaxis()->CenterTitle(kTRUE);
    h_distance->GetYaxis()->SetTitleOffset(1.5);
    h_distance->Draw("HIST");
    h_sftm_distance->Draw("HIST SAME");
    h_ml_distance->Draw("HIST SAME");
    h_wpcluster_distance->Draw("HIST SAME");
    h_cdcluster_distance->Draw("HIST SAME");

    TLegend* leg_distance = new TLegend(0.55, 0.65, 0.85, 0.85);
    leg_distance->AddEntry(h_distance, "Joint #chi^2", "l");
    leg_distance->AddEntry(h_sftm_distance, "SFTM", "l");
    leg_distance->AddEntry(h_ml_distance, "ML", "l");
    leg_distance->AddEntry(h_wpcluster_distance, "WP cluster", "l");
    leg_distance->AddEntry(h_cdcluster_distance, "CD cluster", "l");
    leg_distance->Draw();

    TLatex* tex1_distance = new TLatex();
    tex1_distance->SetNDC();
    tex1_distance->SetTextFont(72);
    tex1_distance->SetTextSize(0.035);
    tex1_distance->SetTextAlign(13); // left-bottom alignment
    tex1_distance->DrawLatex(x_text, y_text, Form("RUN %d-%d", min_run, max_run));

    TLatex* tex2_distance = new TLatex();
    tex2_distance->SetNDC();
    tex2_distance->SetTextFont(42);
    tex2_distance->SetTextSize(0.025);
    tex2_distance->SetTextAlign(13); // left-bottom alignment
    tex2_distance->DrawLatex(x_text, y_text - 0.04, Form("Entries %lld", static_cast<long long>(h_angle->GetEntries())));

    c_distance->SetTickx();
    c_distance->SetTicky();
    c_distance->Update();

    TCanvas* c_angle_clippingness = new TCanvas("c_angle_clippingness", "c_angle_clippingness", 1000, 1000);
    c_angle_clippingness->cd();
    h_angle_clippingness->SetStats(0);
    h_angle_clippingness->GetXaxis()->SetTitle("#alpha (deg)");
    h_angle_clippingness->GetXaxis()->CenterTitle(kTRUE);
    h_angle_clippingness->GetYaxis()->SetTitle("L (m)");
    h_angle_clippingness->GetYaxis()->CenterTitle(kTRUE);
    h_angle_clippingness->GetYaxis()->SetTitleOffset(1.25);
    h_angle_clippingness->Draw("COLZ");
    c_angle_clippingness->Update();

    TCanvas* c_distance_clippingness = new TCanvas("c_distance_clippingness", "c_distance_clippingness", 1000, 1000);
    c_distance_clippingness->cd();
    h_distance_clippingness->SetStats(0);
    h_distance_clippingness->GetXaxis()->SetTitle("d_{mid} (m)");
    h_distance_clippingness->GetXaxis()->CenterTitle(kTRUE);
    h_distance_clippingness->GetYaxis()->SetTitle("L (m)");
    h_distance_clippingness->GetYaxis()->CenterTitle(kTRUE);
    h_distance_clippingness->GetYaxis()->SetTitleOffset(1.25);
    h_distance_clippingness->Draw("COLZ");
    c_distance_clippingness->Update();

    TCanvas* c_angle_distance = new TCanvas("c_angle_distance", "c_angle_distance", 1000, 1000);
    c_angle_distance->cd();
    h_angle_distance->SetStats(0);
    h_angle_distance->GetXaxis()->SetTitle("#alpha (deg)");
    h_angle_distance->GetXaxis()->CenterTitle(kTRUE);
    h_angle_distance->GetYaxis()->SetTitle("d_{mid} (m)");
    h_angle_distance->GetYaxis()->CenterTitle(kTRUE);
    h_angle_distance->GetYaxis()->SetTitleOffset(1.25);
    h_angle_distance->Draw("COLZ");
    c_angle_distance->Update();

    TCanvas* c_angle_chi2cdwp = new TCanvas("c_angle_chi2cdwp", "c_angle_chi2cdwp", 1000, 1000);
    c_angle_chi2cdwp->cd();
    h_angle_chi2cdwp->Draw("COLZ");
    c_angle_chi2cdwp->Update();

    TCanvas* c_distance_chi2cdwp = new TCanvas("c_distance_chi2cdwp", "c_distance_chi2cdwp", 1000, 1000);
    c_distance_chi2cdwp->cd();
    h_distance_chi2cdwp->Draw("COLZ");
    c_distance_chi2cdwp->Update();

    std::nth_element(angles.begin(), angles.begin() + angles.size() * 682 / 1000, angles.end());
    std::nth_element(distances.begin(), distances.begin() + distances.size() * 682 / 1000, distances.end());
    std::cout << "68.2% angle: " << angles[angles.size() * 682 / 1000] << ", size: " << angles.size() << '\n';
    std::cout << "68.2% distance: " << distances[distances.size() * 682 / 1000] << ", size: " << distances.size() << '\n';


    double probs[3] = {0.50, 0.682, 0.90};   // median, 68.2%, 90%
    double q_joint[3], q_sftm[3], q_ml[3], q_wp[3], q_cd[3];

    h_angle->GetQuantiles(3, q_joint, probs);
    h_sftm_angle->GetQuantiles(3, q_sftm, probs);
    h_ml_angle->GetQuantiles(3, q_ml, probs);
    h_wpcluster_angle->GetQuantiles(3, q_wp, probs);
    h_cdcluster_angle->GetQuantiles(3, q_cd, probs);

    std::cout << "\n=== ANGLE QUANTILES ===\n";
    std::cout << "Joint  : 50%=" << q_joint[0] << "  68.2%=" << q_joint[1] << "  90%=" << q_joint[2] << "\n";
    std::cout << "SFTM   : 50%=" << q_sftm[0] << "  68.2%=" << q_sftm[1] << "  90%=" << q_sftm[2] << "\n";
    std::cout << "ML     : 50%=" << q_ml[0]  << "  68.2%=" << q_ml[1]  << "  90%=" << q_ml[2]  << "\n";
    std::cout << "WP cls.: 50%=" << q_wp[0]  << "  68.2%=" << q_wp[1]  << "  90%=" << q_wp[2]  << "\n";
    std::cout << "CD cls.: 50%=" << q_cd[0]  << "  68.2%=" << q_cd[1]  << "  90%=" << q_cd[2]  << "\n";

    h_distance->GetQuantiles(3, q_joint, probs);
    h_sftm_distance->GetQuantiles(3, q_sftm, probs);
    h_ml_distance->GetQuantiles(3, q_ml, probs);
    h_wpcluster_distance->GetQuantiles(3, q_wp, probs);
    h_cdcluster_distance->GetQuantiles(3, q_cd, probs);

    std::cout << "\n=== DISTANCE QUANTILES ===\n";
    std::cout << "Joint  : 50%=" << q_joint[0] << "  68.2%=" << q_joint[1] << "  90%=" << q_joint[2] << "\n";
    std::cout << "SFTM   : 50%=" << q_sftm[0] << "  68.2%=" << q_sftm[1] << "  90%=" << q_sftm[2] << "\n";
    std::cout << "ML     : 50%=" << q_ml[0]  << "  68.2%=" << q_ml[1]  << "  90%=" << q_ml[2]  << "\n";
    std::cout << "WP cls.: 50%=" << q_wp[0]  << "  68.2%=" << q_wp[1]  << "  90%=" << q_wp[2]  << "\n";
    std::cout << "CD cls.: 50%=" << q_cd[0]  << "  68.2%=" << q_cd[1]  << "  90%=" << q_cd[2]  << "\n";

    TH1D* h_sim_tt_zenith = new TH1D("h_sim_tt_zenith", "h_sim_tt_zenith", 100, 0.0, 1.0);
    TH1D* h_sim_tt_azimuth = new TH1D("h_sim_tt_azimuth", "h_sim_tt_azimuth", 100, 0.0, 360.0);

    // TFile* simfile = TFile::Open("/home/traymond/Documents/test/reconstruction_summary/muons.TT_simulation.root", "READ");
    // TTree* simmu = simfile->Get<TTree>("muon");
    // int multiplicity;
    // std::vector<int> *status = nullptr, *pdg = nullptr;
    // std::vector<double> *mass = nullptr, *px = nullptr, *py = nullptr, *pz = nullptr, *x = nullptr, *y = nullptr, *z = nullptr;
    // simmu->SetBranchAddress("multiplicity", &multiplicity);
    // simmu->SetBranchAddress("status", &status);
    // simmu->SetBranchAddress("pdg", &pdg);
    // simmu->SetBranchAddress("px", &px);
    // simmu->SetBranchAddress("py", &py);
    // simmu->SetBranchAddress("pz", &pz);
    // simmu->SetBranchAddress("mass", &mass);
    // simmu->SetBranchAddress("x", &x);
    // simmu->SetBranchAddress("y", &y);
    // simmu->SetBranchAddress("z", &z);

    // double xmin = -23500, xmax = 23500;
    // double ymin = -10000, ymax = 10000;
    // double z1 = 24500;
    // double z2 = 27500;

    // for (int k = 0; k < simmu->GetEntries(); ++k) {
    //     simmu->GetEntry(k);
    //     if (multiplicity != 1) continue;

    //     TVector3 pos((*x)[0], (*y)[0], (*z)[0]);
    //     TVector3 dir((*px)[0], (*py)[0], (*pz)[0]);
    //     dir = dir.Unit();

    //     if (std::abs(dir.Z()) < 1e-9) continue;

    //     auto intersectsPlane = [&](double zplane, double &ix, double &iy) {
    //         double t = (zplane - pos.Z()) / dir.Z();
    //         if (t <= 0) return false;  // behind the muon
    //         ix = pos.X() + t * dir.X();
    //         iy = pos.Y() + t * dir.Y();
    //         return (ix >= xmin && ix <= xmax && iy >= ymin && iy <= ymax);
    //     };

    //     double x1, y1, x2, y2;
    //     bool hit1 = intersectsPlane(z1, x1, y1);
    //     bool hit2 = intersectsPlane(z2, x2, y2);

    //     if (!hit1 || !hit2) continue;

    //     h_sim_tt_zenith->Fill(-dir.CosTheta());
    //     double phi_deg = dir.Phi() * 180.0 / M_PI;
    //     if (phi_deg < 0) phi_deg += 360.0;
    //     h_sim_tt_azimuth->Fill(phi_deg);
    // }

    TFile* simfile = TFile::Open("/home/traymond/Documents/test/reconstruction_summary/rec.TT_simreco.user.root", "READ");
    TTree* ttrecomu = simfile->Get<TTree>("TT");
    ttrecomu->SetBranchAddress("NTotPoints", &NTotPoints);
    ttrecomu->SetBranchAddress("PointX", &PointX);
    ttrecomu->SetBranchAddress("PointY", &PointY);
    ttrecomu->SetBranchAddress("PointZ", &PointZ);
    ttrecomu->SetBranchAddress("NTracks", &NTracks);
    ttrecomu->SetBranchAddress("NPoints", NPoints);
    ttrecomu->SetBranchAddress("start_TS", &start_TS);
    ttrecomu->SetBranchAddress("Coeff0", &Coeff0);
    ttrecomu->SetBranchAddress("Coeff1", &Coeff1);
    ttrecomu->SetBranchAddress("Coeff2", &Coeff2);
    ttrecomu->SetBranchAddress("Coeff3", &Coeff3);
    ttrecomu->SetBranchAddress("Coeff4", &Coeff4);
    ttrecomu->SetBranchAddress("Coeff5", &Coeff5);
    ttrecomu->SetBranchAddress("Chi2", &Chi2);

    for (int k = 0; k < ttrecomu->GetEntries(); ++k) {
        ttrecomu->GetEntry(k);
        if (NTracks != 1) continue;
        TVector3 pos(Coeff0[0], Coeff1[0], Coeff2[0] + cdwp_tt_coordinate_offset);
        TVector3 dir = TVector3(Coeff3[0], Coeff4[0], Coeff5[0]).Unit();
        double distance = dir.Cross(-pos).Mag();
        if (distance >= 17700.0) continue;
        h_sim_tt_zenith->Fill(-dir.CosTheta());
        h_sim_tt_azimuth->Fill( 180.0 / M_PI * (dir.Phi() > 0.0 ? dir.Phi() : dir.Phi() + 2.0 * M_PI) );
    }

    h_sim_tt_zenith->Scale(1.0 / h_sim_tt_zenith->Integral());
    h_sim_tt_azimuth->Scale(1.0 / h_sim_tt_azimuth->Integral());

    TH1D* h_cdwp_zenith = new TH1D("h_cdwp_zenith", "h_cdwp_zenith", 100, 0.0, 1.0);
    TH1D* h_cdwp_azimuth = new TH1D("h_cdwp_azimuth", "h_cdwp_azimuth", 100, 0.0, 360.0);
    TH1D* h_cdwp_zenith_good = new TH1D("h_cdwp_zenith_good", "h_cdwp_zenith_good", 100, 0.0, 1.0);
    TH1D* h_cdwp_azimuth_good = new TH1D("h_cdwp_azimuth_good", "h_cdwp_azimuth_good", 100, 0.0, 360.0);
    TH1D* h_tt_zenith = new TH1D("h_tt_zenith", "h_tt_zenith", 100, 0.0, 1.0);
    TH1D* h_tt_azimuth = new TH1D("h_tt_azimuth", "h_tt_azimuth", 100, 0.0, 360.0);

    for (std::size_t k = 0ul; k < zenith_cdwp_v.size(); ++k) {
        h_cdwp_zenith->Fill(zenith_cdwp_v[k]);
        h_cdwp_azimuth->Fill(azimuth_cdwp_v[k]);
        if (angles[k] < 5.0 && distances[k] < 2.0) {
            h_cdwp_zenith_good->Fill(zenith_cdwp_v[k]);
            h_cdwp_azimuth_good->Fill(azimuth_cdwp_v[k]);
        }
        h_tt_zenith->Fill(zenith_tt_v[k]);
        h_tt_azimuth->Fill(azimuth_tt_v[k]);
    }
    h_cdwp_zenith->Scale(1.0 / h_cdwp_zenith->Integral());
    h_cdwp_azimuth->Scale(1.0 / h_cdwp_azimuth->Integral());
    h_cdwp_zenith_good->Scale(1.0 / h_cdwp_zenith_good->Integral());
    h_cdwp_azimuth_good->Scale(1.0 / h_cdwp_azimuth_good->Integral());
    h_tt_zenith->Scale(1.0 / h_tt_zenith->Integral());
    h_tt_azimuth->Scale(1.0 / h_tt_azimuth->Integral());

    // ============================================================================================
    // Zenith and Azimuth - All correlated "single" muon reconstruction
    // ============================================================================================

    TCanvas* c_sim_tt_zenith = new TCanvas("c_sim_tt_zenith", "c_sim_tt_zenith", 1000, 1000);
    c_sim_tt_zenith->cd();

    h_sim_tt_zenith->SetStats(0);
    h_sim_tt_zenith->SetMinimum(0.0001);
    h_sim_tt_zenith->SetMaximum(std::max({h_sim_tt_zenith->GetMaximum(), h_cdwp_zenith->GetMaximum(), h_tt_zenith->GetMaximum()}));
    h_sim_tt_zenith->GetXaxis()->SetTitle("cos(#theta_{d})");
    h_sim_tt_zenith->GetXaxis()->CenterTitle(kTRUE);
    h_sim_tt_zenith->GetYaxis()->SetTitle("Entries");
    h_sim_tt_zenith->GetYaxis()->CenterTitle(kTRUE);
    h_sim_tt_zenith->GetYaxis()->SetTitleOffset(1.5);
    h_sim_tt_zenith->SetLineWidth(3);
    h_sim_tt_zenith->SetLineColor(kBlack);
    h_sim_tt_zenith->SetLineStyle(kSolid);
    h_sim_tt_zenith->Draw("HIST");
    h_cdwp_zenith->SetLineColor(kRed+1);
    h_cdwp_zenith->SetLineWidth(1);
    h_cdwp_zenith->SetMarkerColor(kRed+1);
    h_cdwp_zenith->SetMarkerSize(1.25);
    h_cdwp_zenith->SetMarkerStyle(kFullCircle);
    h_cdwp_zenith->Draw("SAME");
    h_tt_zenith->SetLineColor(kBlue+1);
    h_tt_zenith->SetLineWidth(1);
    h_tt_zenith->SetMarkerColor(kBlue+1);
    h_tt_zenith->SetMarkerSize(1.25);
    h_tt_zenith->SetMarkerStyle(kFullCircle);
    h_tt_zenith->Draw("SAME");

    TLegend* leg_sim_tt_zenith = new TLegend(0.15, 0.65, 0.45, 0.85);
    leg_sim_tt_zenith->AddEntry(h_sim_tt_zenith, "TT simu reco", "l");
    leg_sim_tt_zenith->AddEntry(h_cdwp_zenith, "CD/WP data reco", "l");
    leg_sim_tt_zenith->AddEntry(h_tt_zenith, "TT data reco", "l");
    leg_sim_tt_zenith->Draw();

    c_sim_tt_zenith->SetTickx();
    c_sim_tt_zenith->SetTicky();
    c_sim_tt_zenith->Update();

    TCanvas* c_sim_tt_azimuth = new TCanvas("c_sim_tt_azimuth", "c_sim_tt_azimuth", 1000, 1000);
    c_sim_tt_azimuth->cd();

    h_sim_tt_azimuth->SetStats(0);
    h_sim_tt_azimuth->SetMinimum(0.0001);
    h_sim_tt_azimuth->SetMaximum(std::max({h_sim_tt_azimuth->GetMaximum(), h_cdwp_azimuth->GetMaximum(), h_tt_azimuth->GetMaximum()}));
    h_sim_tt_azimuth->GetXaxis()->SetTitle("#phi_{d}");
    h_sim_tt_azimuth->GetXaxis()->CenterTitle(kTRUE);
    h_sim_tt_azimuth->GetYaxis()->SetTitle("Entries");
    h_sim_tt_azimuth->GetYaxis()->CenterTitle(kTRUE);
    h_sim_tt_azimuth->GetYaxis()->SetTitleOffset(1.5);
    h_sim_tt_azimuth->SetLineWidth(3);
    h_sim_tt_azimuth->SetLineColor(kBlack);
    h_sim_tt_azimuth->SetLineStyle(kSolid);
    h_sim_tt_azimuth->Draw("HIST");
    h_cdwp_azimuth->SetLineColor(kRed+1);
    h_cdwp_azimuth->SetLineWidth(1);
    h_cdwp_azimuth->SetMarkerColor(kRed+1);
    h_cdwp_azimuth->SetMarkerSize(1.25);
    h_cdwp_azimuth->SetMarkerStyle(kFullCircle);
    h_cdwp_azimuth->Draw("SAME");
    h_tt_azimuth->SetLineColor(kBlue+1);
    h_tt_azimuth->SetLineWidth(1);
    h_tt_azimuth->SetMarkerColor(kBlue+1);
    h_tt_azimuth->SetMarkerSize(1.25);
    h_tt_azimuth->SetMarkerStyle(kFullCircle);
    h_tt_azimuth->Draw("SAME");

    TLegend* leg_sim_tt_azimuth = new TLegend(0.35, 0.15, 0.65, 0.35);
    leg_sim_tt_azimuth->AddEntry(h_sim_tt_azimuth, "TT simu reco", "l");
    leg_sim_tt_azimuth->AddEntry(h_cdwp_azimuth, "CD/WP data reco", "l");
    leg_sim_tt_azimuth->AddEntry(h_tt_azimuth, "TT data reco", "l");
    leg_sim_tt_azimuth->Draw();

    c_sim_tt_azimuth->SetTickx();
    c_sim_tt_azimuth->SetTicky();
    c_sim_tt_azimuth->Update();

    // ============================================================================================
    // Zenith and Azimuth - Only "good" correlated "single" muon reconstruction
    // ============================================================================================

    TCanvas* c_sim_tt_zenith_good = new TCanvas("c_sim_tt_zenith_good", "c_sim_tt_zenith_good", 1000, 1000);
    c_sim_tt_zenith_good->cd();

    h_sim_tt_zenith->SetStats(0);
    h_sim_tt_zenith->SetMinimum(0.0001);
    h_sim_tt_zenith->SetMaximum(std::max({h_sim_tt_zenith->GetMaximum(), h_cdwp_zenith_good->GetMaximum(), h_tt_zenith->GetMaximum()}));
    h_sim_tt_zenith->GetXaxis()->SetTitle("cos(#theta_{d})");
    h_sim_tt_zenith->GetXaxis()->CenterTitle(kTRUE);
    h_sim_tt_zenith->GetYaxis()->SetTitle("Entries");
    h_sim_tt_zenith->GetYaxis()->CenterTitle(kTRUE);
    h_sim_tt_zenith->GetYaxis()->SetTitleOffset(1.5);
    h_sim_tt_zenith->SetLineWidth(3);
    h_sim_tt_zenith->SetLineColor(kBlack);
    h_sim_tt_zenith->SetLineStyle(kSolid);
    h_sim_tt_zenith->Draw("HIST");
    h_cdwp_zenith_good->SetLineColor(kRed+1);
    h_cdwp_zenith_good->SetLineWidth(1);
    h_cdwp_zenith_good->SetMarkerColor(kRed+1);
    h_cdwp_zenith_good->SetMarkerSize(1.25);
    h_cdwp_zenith_good->SetMarkerStyle(kFullCircle);
    h_cdwp_zenith_good->Draw("SAME");
    h_tt_zenith->SetLineColor(kBlue+1);
    h_tt_zenith->SetLineWidth(1);
    h_tt_zenith->SetMarkerColor(kBlue+1);
    h_tt_zenith->SetMarkerSize(1.25);
    h_tt_zenith->SetMarkerStyle(kFullCircle);
    h_tt_zenith->Draw("SAME");

    TLegend* leg_sim_tt_zenith_good = new TLegend(0.15, 0.65, 0.45, 0.85);
    leg_sim_tt_zenith_good->AddEntry(h_sim_tt_zenith, "TT simu reco", "l");
    leg_sim_tt_zenith_good->AddEntry(h_cdwp_zenith_good, "CD/WP data reco", "l");
    leg_sim_tt_zenith_good->AddEntry(h_tt_zenith, "TT data reco", "l");
    leg_sim_tt_zenith_good->Draw();

    c_sim_tt_zenith_good->SetTickx();
    c_sim_tt_zenith_good->SetTicky();
    c_sim_tt_zenith_good->Update();

    TCanvas* c_sim_tt_azimuth_good = new TCanvas("c_sim_tt_azimuth_good", "c_sim_tt_azimuth_good", 1000, 1000);
    c_sim_tt_azimuth_good->cd();

    h_sim_tt_azimuth->SetStats(0);
    h_sim_tt_azimuth->SetMinimum(0.0001);
    h_sim_tt_azimuth->SetMaximum(std::max({h_sim_tt_azimuth->GetMaximum(), h_cdwp_azimuth_good->GetMaximum(), h_tt_azimuth->GetMaximum()}));
    h_sim_tt_azimuth->GetXaxis()->SetTitle("#phi_{d}");
    h_sim_tt_azimuth->GetXaxis()->CenterTitle(kTRUE);
    h_sim_tt_azimuth->GetYaxis()->SetTitle("Entries");
    h_sim_tt_azimuth->GetYaxis()->CenterTitle(kTRUE);
    h_sim_tt_azimuth->GetYaxis()->SetTitleOffset(1.5);
    h_sim_tt_azimuth->SetLineWidth(3);
    h_sim_tt_azimuth->SetLineColor(kBlack);
    h_sim_tt_azimuth->SetLineStyle(kSolid);
    h_sim_tt_azimuth->Draw("HIST");
    h_cdwp_azimuth_good->SetLineColor(kRed+1);
    h_cdwp_azimuth_good->SetLineWidth(1);
    h_cdwp_azimuth_good->SetMarkerColor(kRed+1);
    h_cdwp_azimuth_good->SetMarkerSize(1.25);
    h_cdwp_azimuth_good->SetMarkerStyle(kFullCircle);
    h_cdwp_azimuth_good->Draw("SAME");
    h_tt_azimuth->SetLineColor(kBlue+1);
    h_tt_azimuth->SetLineWidth(1);
    h_tt_azimuth->SetMarkerColor(kBlue+1);
    h_tt_azimuth->SetMarkerSize(1.25);
    h_tt_azimuth->SetMarkerStyle(kFullCircle);
    h_tt_azimuth->Draw("SAME");

    TLegend* leg_sim_tt_azimuth_good = new TLegend(0.35, 0.15, 0.65, 0.35);
    leg_sim_tt_azimuth_good->AddEntry(h_sim_tt_azimuth, "TT simu reco", "l");
    leg_sim_tt_azimuth_good->AddEntry(h_cdwp_azimuth_good, "CD/WP data reco", "l");
    leg_sim_tt_azimuth_good->AddEntry(h_tt_azimuth, "TT data reco", "l");
    leg_sim_tt_azimuth_good->Draw();

    c_sim_tt_azimuth_good->SetTickx();
    c_sim_tt_azimuth_good->SetTicky();
    c_sim_tt_azimuth_good->Update();
}