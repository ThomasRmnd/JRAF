#include <cstdint>
#include <iostream>

#include <TFile.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TVector3.h>

struct track {

    int run_id;
    TTimeStamp ts;
    double totq_cd;
    double totq_wp;
    double quality;

    TVector3 ipos;
    TVector3 fpos;

};

bool operator<(const track& a, const track& b) {
    return a.ts < b.ts;
}

double compute_angle_between_track(const track& a, const track& b) {
    TVector3 dir_a = (a.fpos - a.ipos).Unit();
    TVector3 dir_b = (b.fpos - b.ipos).Unit();
    return dir_a.Angle(dir_b) * 180.0 / M_PI;
}

double compute_distance_between_track(const track& a, const track& b) {
    TVector3 mpos_a = (a.ipos + a.fpos) * 0.5;
    TVector3 mpos_b = (b.ipos + b.fpos) * 0.5;
    return (mpos_a - mpos_b).Mag() / 1000.0;
}

double compute_clippingness(const track& a) {
    TVector3 dir = (a.fpos - a.ipos).Unit();
    return dir.Cross(-a.ipos).Mag() / dir.Mag() / 1000.0;
}

double get_quantile(std::vector<double>::const_iterator first, std::vector<double>::const_iterator last, double quantile) {
    if (first == last) return 0.0;
    std::vector<double> v(first, last);
    std::vector<double>::iterator position = v.begin() + quantile * v.size();
    std::nth_element(v.begin(), position, v.end());
    return *position;
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

std::set<track> open_amber_v5_5_user_chain(const char* path) {
    TChain* chain = new TChain("MuonReco");
    chain->Add(path);
    std::set<track> tracks;

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
    std::size_t nstoppins = 0ul;
    std::cout << "Info: Found " << nentries << " entries in Amber_v5.5 files\n";
    for (long k = 0l; k < nentries; ++k) {
        chain->GetEntry(k);
        if (muonType != 0) continue; // SELECTION! only single
        TVector3 fpos(xout, yout, zout);
        if (fpos.Mag() > 40000.0) {
            ++nstoppins;
        }
        tracks.insert(track{
            .run_id = runID,
            .ts = TTimeStamp(fSec, fNanoSec),
            .totq_cd = 0.0,
            .totq_wp = charge,
            .quality = 0.0,
            .ipos = TVector3(xin, yin, zin),
            .fpos = fpos
        });
    }
    std::cout << "Info: Number of stopping tracks for Amber: " << nstoppins << '\n';
    return tracks;
}

std::set<track> open_edwin_user_chain(const char* path) {
    TChain* chain = new TChain("Single_Reco");
    chain->Add(path);
    std::set<track> tracks;

    int run_number;
    int cd_file;
    int cd_time_s;
    int cd_time_ns;
    float cd_totalPE;
    float enterX, enterY, enterZ;
    float exitX, exitY, exitZ;
    int wp_time_s;
    int wp_time_ns;
    float wp_totalPE;
    chain->SetBranchAddress("run_number", &run_number);
    chain->SetBranchAddress("cd_file", &cd_file);
    chain->SetBranchAddress("cd_time_s", &cd_time_s);
    chain->SetBranchAddress("cd_time_ns", &cd_time_ns);
    chain->SetBranchAddress("cd_totalPE", &cd_totalPE);
    chain->SetBranchAddress("enterX", &enterX);
    chain->SetBranchAddress("enterY", &enterY);
    chain->SetBranchAddress("enterZ", &enterZ);
    chain->SetBranchAddress("exitX", &exitX);
    chain->SetBranchAddress("exitY", &exitY);
    chain->SetBranchAddress("exitZ", &exitZ);
    chain->SetBranchAddress("wp_time_s", &wp_time_s);
    chain->SetBranchAddress("wp_time_ns", &wp_time_ns);
    chain->SetBranchAddress("wp_totalPE", &wp_totalPE);

    long nentries = chain->GetEntries();
    std::size_t nstoppins = 0ul;
    std::cout << "Info: Found " << nentries << " entries in EDWIN files\n";
    for (long k = 0l; k < nentries; ++k) {
        chain->GetEntry(k);
        TVector3 fpos(exitX, exitY, exitZ);
        if (fpos.Mag() > 40000.0) {
            ++nstoppins;
        }
        tracks.insert(track{
            .run_id = run_number,
            .ts = TTimeStamp(cd_time_s, cd_time_ns),
            .totq_cd = cd_totalPE,
            .totq_wp = wp_totalPE,
            .quality = 0.0,
            .ipos = TVector3(enterX, enterY, enterZ),
            .fpos = fpos
        });
    }
    std::cout << "Info: Number of stopping tracks for Edwin: " << nstoppins << '\n';
    return tracks;
}

std::set<track> open_cdwpttchi2_user_chain(const char* path) {
    TChain* chain = new TChain("muons");
    chain->Add(path);
    std::set<track> tracks;

    int run_id;
    time_t sec;
    int nsec;
    double totq_cd;
    double totq_wp;
    double chi2;
    double iposx, iposy, iposz;
    double fposx, fposy, fposz;

    chain->SetBranchAddress("run_id", &run_id);
    chain->SetBranchAddress("sec", &sec);
    chain->SetBranchAddress("nsec", &nsec);
    chain->SetBranchAddress("totq_cd", &totq_cd);
    chain->SetBranchAddress("totq_wp", &totq_wp);
    chain->SetBranchAddress("chi2", &chi2);
    chain->SetBranchAddress("iposx", &iposx);
    chain->SetBranchAddress("iposy", &iposy);
    chain->SetBranchAddress("iposz", &iposz);
    chain->SetBranchAddress("fposx", &fposx);
    chain->SetBranchAddress("fposy", &fposy);
    chain->SetBranchAddress("fposz", &fposz);

    long nentries = chain->GetEntries();
    std::cout << "Info: Found " << nentries << " entries in CdWpTtChi2 files\n";
    for (long k = 0l; k < nentries; ++k) {
        chain->GetEntry(k);
        tracks.insert(track{
            .run_id = run_id,
            .ts = TTimeStamp(sec, nsec),
            .totq_cd = totq_cd,
            .totq_wp = totq_wp,
            .quality = chi2,
            .ipos = TVector3(iposx, iposy, iposz),
            .fpos = TVector3(fposx, fposy, fposz)
        });
    }
    return tracks;
}

std::map<std::string, std::set<track>> open_joint_reco_user_chain(const char* path) {
    TChain* chain = new TChain("muons");
    chain->Add(path);
    std::map<std::string, std::set<track>> tracks;

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
        bool has_tt_info = false;
        bool is_in_acrylic = false;
        int ntracks_cdclassify = 0;
        bool stopping_cdclassify = false;
        int ntracks_wpclassify = 0;
        bool stopping_wpclassify = false;
        for (std::size_t i = 0ul; i < method->size(); ++i) {
            if ((*method)[i] == "Tt") {
                has_tt_info = true;
                TVector3 ipos((*iposx)[i], (*iposy)[i], (*iposz)[i]);
                TVector3 fpos((*fposx)[i], (*fposy)[i], (*fposz)[i]);
                TVector3 dir = (fpos - ipos).Unit();
                if (dir.Cross(-ipos).Mag() < 17700.0) {
                    is_in_acrylic = true;
                }
            }
            if ((*method)[i] == "CdClassify") {
                ++ntracks_cdclassify;
                TVector3 ipos((*iposx)[i], (*iposy)[i], (*iposz)[i]);
                TVector3 fpos((*fposx)[i], (*fposy)[i], (*fposz)[i]);
                if (fpos.Mag() > 40000.0) {
                    stopping_cdclassify = true;
                }
            }
            if ((*method)[i] == "WpBasic") {
                ++ntracks_wpclassify;

                TVector3 ipos((*iposx)[i], (*iposy)[i], (*iposz)[i]);
                TVector3 fpos((*fposx)[i], (*fposy)[i], (*fposz)[i]);
                if (fpos.Mag() > 40000.0) {
                    stopping_wpclassify = true;
                }
            }
        }
        if (!has_tt_info) continue;
        // if (!is_in_acrylic) continue; // SELECTION!
        if (ntracks_cdclassify != 1) continue; // SELECTION! || stopping_cdclassify
        if (ntracks_wpclassify != 1) continue; // SELECTION! || stopping_wpclassify
        for (std::size_t i = 0ul; i < method->size(); ++i) {
            if ((*method)[i] == "CdWpTtChi2") continue;
            tracks[(*method)[i]].insert(track{
                .run_id = run_id,
                .ts = TTimeStamp(sec, nsec),
                .totq_cd = totq_cd,
                .totq_wp = totq_wp,
                .quality = (*quality)[i],
                .ipos = TVector3((*iposx)[i], (*iposy)[i], (*iposz)[i]),
                .fpos = TVector3((*fposx)[i], (*fposy)[i], (*fposz)[i])
            });
        }
    }
    return tracks;
}

struct MuonPerformance {
    double angle;
    double distance;
    double clippingness;
    double quality;
    double tt_quality;
};

std::vector<double> extract_angles_from_performances(const std::vector<MuonPerformance>& perf) {
    std::vector<double> angles;
    for (const MuonPerformance& mp : perf) {
        angles.push_back(mp.angle);
    }
    return angles;
}

std::vector<double> extract_distances_from_performances(const std::vector<MuonPerformance>& perf) {
    std::vector<double> distances;
    for (const MuonPerformance& mp : perf) {
        distances.push_back(mp.distance);
    }
    return distances;
}

std::map<std::string, std::vector<MuonPerformance>> compute_correlations(std::map<std::string, std::set<track>>& tracks) {
    if (tracks.find("Tt") == tracks.end()) {
        std::cerr << "Error: Tt tracks not found in map.\n";
        return {};
    }
    const std::set<track>& tt_tracks = tracks["Tt"];

    std::map<std::string, std::vector<MuonPerformance>> performances;

    for (const auto& [method, track_set] : tracks) {
        if (method == "Tt") continue;

        std::cout << "\n--- Correlating " << method << " with Tt ---" << std::endl;

        for (const track& trk : track_set) {
            TTimeStamp lower_bound_ts(trk.ts.GetSec(), trk.ts.GetNanoSec() - 1000);
            TTimeStamp upper_bound_ts(trk.ts.GetSec(), trk.ts.GetNanoSec() + 1000);
            
            std::set<track>::const_iterator it_tt = tt_tracks.lower_bound({0, lower_bound_ts, 0, 0, {}, {}});

            while (it_tt != track_set.end() && lower_bound_ts <= it_tt->ts && it_tt->ts <= upper_bound_ts) {
                performances[method].push_back(MuonPerformance{
                    .angle = compute_angle_between_track(trk, *it_tt),
                    .distance = compute_distance_between_track(trk, *it_tt),
                    .clippingness = compute_clippingness(*it_tt),
                    .quality = trk.quality,
                    .tt_quality = it_tt->quality
                });
                ++it_tt;
            }
        }
    }
    return performances;
}

std::map<std::string, std::vector<MuonPerformance>> compute_global_correlations(std::map<std::string, std::set<track>>& tracks) {
    if (tracks.find("Tt") == tracks.end()) {
        std::cerr << "Error: Tt tracks not found in map.\n";
        return {};
    }
    const std::set<track>& tt_tracks = tracks["Tt"];

    std::map<std::string, std::vector<MuonPerformance>> performances;

    for (const track& tt_muon : tt_tracks) {
        std::map<std::string, track> coincident_map;
        bool all_found = true;

        for (const auto& [method, track_set] : tracks) {
            if (method == "Tt") continue;
            bool found_in_method = false;
            TTimeStamp lower_bound_ts(tt_muon.ts.GetSec(), tt_muon.ts.GetNanoSec() - 1000);
            TTimeStamp upper_bound_ts(tt_muon.ts.GetSec(), tt_muon.ts.GetNanoSec() + 1000);
            std::set<track>::const_iterator it = track_set.lower_bound({0, lower_bound_ts, 0, 0, {}, {}});
            
            while (it != track_set.end() && lower_bound_ts <= it->ts && it->ts <= upper_bound_ts) {
                coincident_map[method] = *it;
                found_in_method = true;
                break;
            }

            if (!found_in_method) {
                all_found = false;
                break; 
            }
        }

        if (all_found) {
            for (const auto& [method, muon] : coincident_map) {
                performances[method].push_back(MuonPerformance{
                    .angle = compute_angle_between_track(tt_muon, muon),
                    .distance = compute_distance_between_track(tt_muon, muon),
                    .clippingness = compute_clippingness(tt_muon),
                    .quality = muon.quality,
                    .tt_quality = tt_muon.quality
                });
            }
        }
    }
    return performances;
}

int fast_muon_reconstruction_comparison(const char* path_joint, const char* path_cdwpttchi2, const char* path_amber, const char* path_edwin) {
    std::map<std::string, std::set<track>> tracks = open_joint_reco_user_chain(path_joint);
    tracks["CdWpTtChi2"] = open_cdwpttchi2_user_chain(path_cdwpttchi2);
    tracks["Amber_v5.5"] = open_amber_v5_5_user_chain(path_amber);
    tracks["Edwin"] = open_edwin_user_chain(path_edwin);

    // std::map<std::string, std::vector<MuonPerformance>> performances = compute_correlations(tracks);
    std::map<std::string, std::vector<MuonPerformance>> performances = compute_global_correlations(tracks);
    
    std::map<std::string, std::vector<double>> angles;
    std::map<std::string, std::vector<double>> distances;
    for (const auto& [method, perf] : performances) {
        angles[method] = extract_angles_from_performances(perf);
        distances[method] = extract_distances_from_performances(perf);
    }

    std::map<std::string, TH1D*> method_angle_map;
    std::map<std::string, TH1D*> method_distance_map;

    double xmin_angle = 0.0, xmax_angle = 180.0;
    double xmin_distance = 0.0, xmax_distance = 40.0;
    int nbins_angle = 200, nbins_distance = 200;
    // double xmin_angle = 0.0, xmax_angle = 5.0;
    // double xmin_distance = 0.0, xmax_distance = 2.0;
    // int nbins_angle = 50, nbins_distance = 50;

    method_angle_map["CdWpTtChi2"] = new TH1D("h_angle_cdwpttchi2", "Angle between tracks direction (CdWpTtChi2);#alpha (deg);Entries;", nbins_angle, xmin_angle, xmax_angle);
    method_distance_map["CdWpTtChi2"] = new TH1D("h_distance_cdwpttchi2", "Distance between tracks middle point (CdWpTtChi2);d_{mid} (m);Entries;", nbins_distance, xmin_distance, xmax_distance);
    method_angle_map["CdClassify"] = new TH1D("h_angle_cdclassify", "Angle between tracks direction (CdClassify);#alpha (deg);Entries;", nbins_angle, xmin_angle, xmax_angle);
    method_distance_map["CdClassify"] = new TH1D("h_distance_cdclassify", "Distance between tracks middle point (CdClassify);d_{mid} (m);Entries;", nbins_distance, xmin_distance, xmax_distance);
    method_angle_map["WpBasic"] = new TH1D("h_angle_wpclassify", "Angle between tracks direction (WpClassify);#alpha (deg);Entries;", nbins_angle, xmin_angle, xmax_angle);
    method_distance_map["WpBasic"] = new TH1D("h_distance_wpclassify", "Distance between tracks middle point (WpClassify);d_{mid} (m);Entries;", nbins_distance, xmin_distance, xmax_distance);
    method_angle_map["Amber_v5.5"] = new TH1D("h_angle_amber", "Angle between tracks direction (Amber);#alpha (deg);Entries;", nbins_angle, xmin_angle, xmax_angle);
    method_distance_map["Amber_v5.5"] = new TH1D("h_distance_amber", "Distance between tracks middle point (Amber);d_{mid} (m);Entries;", nbins_distance, xmin_distance, xmax_distance);
    method_angle_map["Edwin"] = new TH1D("h_angle_edwin", "Angle between tracks direction (Edwin);#alpha (deg);Entries;", nbins_angle, xmin_angle, xmax_angle);
    method_distance_map["Edwin"] = new TH1D("h_distance_edwin", "Distance between tracks middle point (Edwin);d_{mid} (m);Entries;", nbins_distance, xmin_distance, xmax_distance);

    for (const auto& [method, perf] : performances) {
        for (const MuonPerformance& mp : perf) {
            method_angle_map[method]->Fill(mp.angle);
            method_distance_map[method]->Fill(mp.distance);
        }
    }

    // std::cout << "68.2% angle: " << get_quantile(angles.begin(), angles.end(), 0.682) << '\n';
    // std::cout << "95.4% angle: " << get_quantile(angles.begin(), angles.end(), 0.954) << '\n';
    // std::cout << "68.2% distance: " << get_quantile(distances.begin(), distances.end(), 0.682) << '\n';
    // std::cout << "95.4% distance: " << get_quantile(distances.begin(), distances.end(), 0.954) << '\n';

    std::map<std::string, double> angle_quantiles;
    for (const auto& [method, ang] : angles) {
        angle_quantiles[method] = get_quantile(ang.begin(), ang.end(), 0.682);
    }
    std::map<std::string, double> distance_quantiles;
    for (const auto& [method, dist] : distances) {
        distance_quantiles[method] = get_quantile(dist.begin(), dist.end(), 0.682);
    }


    plot_metrics(method_angle_map, angle_quantiles);
    plot_metrics(method_distance_map, distance_quantiles);

    double r2_min = 0.0 * 0.0, r2_max = 18.0 * 18.0;
    int nbins = 9;
    std::map<std::string, std::vector<std::vector<double>>> method_angle_r2_bin_content;
    std::map<std::string, std::vector<std::vector<double>>> method_distance_r2_bin_content;

    std::map<std::string, TH1D*> method_angle_r2_map;
    method_angle_r2_map["CdWpTtChi2"] = new TH1D("h_angle_r2_cdwpttchi2", "Angle between tracks direction (CdWpTtChi2);L (m); 68% quantile of #alpha (deg);", nbins, r2_min, r2_max);
    method_angle_r2_map["CdClassify"] = new TH1D("h_angle_r2_cdclassify", "Angle between tracks direction (CdClassify);L (m); 68% quantile of #alpha (deg);", nbins, r2_min, r2_max);
    method_angle_r2_map["WpBasic"] = new TH1D("h_angle_r2_wpclassify", "Angle between tracks direction (WpClassify);L (m); 68% quantile of #alpha (deg);", nbins, r2_min, r2_max);
    method_angle_r2_map["Amber_v5.5"] = new TH1D("h_angle_r2_amber", "Angle between tracks direction (Amber);L (m); 68% quantile of #alpha (deg);", nbins, r2_min, r2_max);
    method_angle_r2_map["Edwin"] = new TH1D("h_angle_r2_edwin", "Angle between tracks direction (Edwin);L (m); 68% quantile of #alpha (deg);", nbins, r2_min, r2_max);
    for (auto& [method, h] : method_angle_r2_map) {
        h->GetXaxis()->SetNdivisions(nbins, false);
    }

    std::map<std::string, TH1D*> method_distance_r2_map;
    method_distance_r2_map["CdWpTtChi2"] = new TH1D("h_distance_r2_cdwpttchi2", "Distance between tracks middle point (CdWpTtChi2);L (m); 68% quantile of d_{mid} (m);", nbins, r2_min, r2_max);
    method_distance_r2_map["CdClassify"] = new TH1D("h_distance_r2_cdclassify", "Distance between tracks middle point (CdClassify);L (m); 68% quantile of d_{mid} (m);", nbins, r2_min, r2_max);
    method_distance_r2_map["WpBasic"] = new TH1D("h_distance_r2_wpclassify", "Distance between tracks middle point (WpClassify);L (m); 68% quantile of d_{mid} (m);", nbins, r2_min, r2_max);
    method_distance_r2_map["Amber_v5.5"] = new TH1D("h_distance_r2_amber", "Distance between tracks middle point (Amber);L (m); 68% quantile of d_{mid} (m);", nbins, r2_min, r2_max);
    method_distance_r2_map["Edwin"] = new TH1D("h_distance_r2_edwin", "Distance between tracks middle point (Edwin);L (m); 68% quantile of d_{mid} (m);", nbins, r2_min, r2_max);
    for (auto& [method, h] : method_distance_r2_map) {
        h->GetXaxis()->SetNdivisions(nbins, false);
    }

    for (const auto& [method, perf] : performances) {
        if (method == "Tt") continue;
        method_angle_r2_bin_content[method].resize(nbins);
        method_distance_r2_bin_content[method].resize(nbins);
        for (const MuonPerformance& mp : perf) {
            double r2 = mp.clippingness * mp.clippingness;
            if (r2 < r2_min || r2_max <= r2) continue;
            int bin = std::floor((r2 - r2_min) / (r2_max - r2_min) * nbins);
            method_angle_r2_bin_content[method][bin].push_back(mp.angle);
            method_distance_r2_bin_content[method][bin].push_back(mp.distance);
        }
        for (int i = 0; i < nbins; ++i) {
            method_angle_r2_map[method]->SetBinContent(i + 1, get_quantile(method_angle_r2_bin_content[method][i].begin(), method_angle_r2_bin_content[method][i].end(), 0.682));
            method_distance_r2_map[method]->SetBinContent(i + 1, get_quantile(method_distance_r2_bin_content[method][i].begin(), method_distance_r2_bin_content[method][i].end(), 0.682));
            method_angle_r2_map[method]->SetBinError(i + 1, 0.0001);
            method_distance_r2_map[method]->SetBinError(i + 1, 0.0001);
            double edge = std::sqrt(i * (r2_max - r2_min) / nbins + r2_min);
            method_angle_r2_map[method]->GetXaxis()->ChangeLabel(i + 1, -1.0, -1.0, -1, -1, -1, Form("%0.1f^{2}", edge));
            method_distance_r2_map[method]->GetXaxis()->ChangeLabel(i + 1, -1.0, -1.0, -1, -1, -1, Form("%0.1f^{2}", edge));
        }
        method_angle_r2_map[method]->GetXaxis()->ChangeLabel(nbins + 1, -1.0, -1.0, -1, -1, -1, Form("%0.1f^{2}", std::sqrt(r2_max)));
        method_distance_r2_map[method]->GetXaxis()->ChangeLabel(nbins + 1, -1.0, -1.0, -1, -1, -1, Form("%0.1f^{2}", std::sqrt(r2_max)));
    }

    std::map<std::string, Color_t> colors = {
        {"CdWpTtChi2", kBlack}, 
        {"CdClassify", kGreen + 2}, 
        {"WpBasic", kViolet}, 
        {"Amber_v5.5", kBlue},
        {"Edwin", kRed},
    };

    TCanvas* c_angle_68p_r2 = new TCanvas("c_angle_68p_r2", "Angle 68% quantile", 1000, 1000);
    c_angle_68p_r2->cd();

    TLegend* leg_angle_68p_r2 = new TLegend(0.15, 0.65, 0.55, 0.85);
    bool is_first_angle = true;

    double max_angle = 0.0;
    for (auto& [method, h] : method_angle_r2_map) {
        if (h->GetMaximum() > max_angle) {
            max_angle = h->GetMaximum();
        }
    }

    for (auto& [method, h] : method_angle_r2_map) {
        h->SetMaximum(0.0);
        h->SetMaximum(10.0);
        h->SetStats(0);
        h->SetMarkerStyle(kFullCircle);
        h->SetMarkerSize(2.0);
        h->SetMarkerColor(colors[method]);
        h->SetLineWidth(2);
        h->SetLineColor(colors[method]);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->GetYaxis()->SetTitleOffset(1.25);
        if (is_first_angle) {
            is_first_angle = false;
            h->Draw("E P");
        }
        else {
            h->Draw("P SAME");
        }
        leg_angle_68p_r2->AddEntry(h, Form("%s: 68%% quantile = %.2f", method.c_str(), angle_quantiles.at(method)), "p");
    }
    leg_angle_68p_r2->SetTextSize(0.02);
    leg_angle_68p_r2->Draw();

    c_angle_68p_r2->SetTickx();
    c_angle_68p_r2->SetTicky();
    c_angle_68p_r2->SetGrid();
    c_angle_68p_r2->Update();

    TCanvas* c_distance_68p_r2 = new TCanvas("c_distance_68p_r2", "Distance 68% quantile", 1000, 1000);
    c_distance_68p_r2->cd();

    double max_distance = 0.0;
    for (auto& [method, h] : method_distance_r2_map) {
        if (h->GetMaximum() > max_distance) {
            max_distance = h->GetMaximum();
        }
    }

    TLegend* leg_distance_68p_r2 = new TLegend(0.15, 0.65, 0.55, 0.85);
    bool is_first_distance = true;

    for (auto& [method, h] : method_distance_r2_map) {
        h->SetMaximum(0.0);
        h->SetMaximum(2.0);
        h->SetStats(0);
        h->SetMarkerStyle(kFullCircle);
        h->SetMarkerSize(2.0);
        h->SetMarkerColor(colors[method]);
        h->SetLineWidth(2);
        h->SetLineColor(colors[method]);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->GetYaxis()->SetTitleOffset(1.25);
        if (is_first_distance) {
            is_first_distance = false;
            h->Draw("E P");
        }
        else {
            h->Draw("P SAME");
        }
        leg_distance_68p_r2->AddEntry(h, Form("%s: 68%% quantile = %.2f", method.c_str(), distance_quantiles.at(method)), "p");
    }
    leg_distance_68p_r2->SetTextSize(0.02);
    leg_distance_68p_r2->Draw();

    c_distance_68p_r2->SetTickx();
    c_distance_68p_r2->SetTicky();
    c_distance_68p_r2->SetGrid();
    c_distance_68p_r2->Update();

    return 0;
}