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

    bool is_single;
    bool is_stopping;

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
        {"Janus", kMagenta},
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
        // if (muonType != 0) continue; // SELECTION! only single
        TVector3 fpos(xout, yout, zout);
        if (fpos.Mag() < 17700.0) {
            ++nstoppins;
        }
        tracks.insert(track{
            .run_id = runID,
            .ts = TTimeStamp(fSec, fNanoSec),
            .totq_cd = 0.0,
            .totq_wp = charge,
            .quality = 0.0,
            .ipos = TVector3(xin, yin, zin),
            .fpos = fpos,
            .is_single = (muonType == 0),
            .is_stopping = false
        });
    }
    std::cout << "Info: Number of stopping tracks for Amber: " << nstoppins << '\n';
    return tracks;
}

std::set<track> open_edwin_user_chain(const char* path) {
    TChain* chain = new TChain("Edwin_Muon");
    chain->Add(path);
    std::set<track> tracks;

    long long cd_time_s;
    long long cd_time_ns;
    int muon_classification; // 0 = single, 1 = double, 2 = triple+
    float Single_enterX;
    float Single_enterY;
    float Single_enterZ;
    float Single_exit_X;
    float Single_exit_Y;
    float Single_exit_Z;
    float Double_enterX_1;
    float Double_enterY_1;
    float Double_enterZ_1;
    float Double_exitX_1;
    float Double_exitY_1;
    float Double_exitZ_1;
    float Double_enterX_2;
    float Double_enterY_2;
    float Double_enterZ_2;
    float Double_exitX_2;
    float Double_exitY_2;
    float Double_exitZ_2;

    chain->SetBranchAddress("cd_time_s", &cd_time_s);
    chain->SetBranchAddress("cd_time_ns", &cd_time_ns);
    chain->SetBranchAddress("muon_classification", &muon_classification);
    chain->SetBranchAddress("Single_enterX", &Single_enterX);
    chain->SetBranchAddress("Single_enterY", &Single_enterY);
    chain->SetBranchAddress("Single_enterZ", &Single_enterZ);
    chain->SetBranchAddress("Single_exit_X", &Single_exit_X);
    chain->SetBranchAddress("Single_exit_Y", &Single_exit_Y);
    chain->SetBranchAddress("Single_exit_Z", &Single_exit_Z);
    chain->SetBranchAddress("Double_enterX_1", &Double_enterX_1);
    chain->SetBranchAddress("Double_enterY_1", &Double_enterY_1);
    chain->SetBranchAddress("Double_enterZ_1", &Double_enterZ_1);
    chain->SetBranchAddress("Double_exitX_1", &Double_exitX_1);
    chain->SetBranchAddress("Double_exitY_1", &Double_exitY_1);
    chain->SetBranchAddress("Double_exitZ_1", &Double_exitZ_1);
    chain->SetBranchAddress("Double_enterX_2", &Double_enterX_2);
    chain->SetBranchAddress("Double_enterY_2", &Double_enterY_2);
    chain->SetBranchAddress("Double_enterZ_2", &Double_enterZ_2);
    chain->SetBranchAddress("Double_exitX_2", &Double_exitX_2);
    chain->SetBranchAddress("Double_exitY_2", &Double_exitY_2);
    chain->SetBranchAddress("Double_exitZ_2", &Double_exitZ_2);

    long nentries = chain->GetEntries();
    std::size_t nstoppins = 0ul;
    std::cout << "Info: Found " << nentries << " entries in EDWIN files\n";
    for (long k = 0l; k < nentries; ++k) {
        chain->GetEntry(k);
        // if (muon_classification != 0) continue; // SELECTION! only single
        tracks.insert(track{
            .run_id = 0,
            .ts = TTimeStamp(static_cast<time_t>(cd_time_s), static_cast<int>(cd_time_ns)),
            .totq_cd = 0.0,
            .totq_wp = 0.0,
            .quality = 0.0,
            .ipos = TVector3(Single_enterX, Single_enterY, Single_enterZ),
            .fpos = TVector3(Single_exit_X, Single_exit_Y, Single_exit_Z),
            .is_single = (muon_classification == 0),
            .is_stopping = false
        });
    }
    return tracks;
}

std::set<track> open_cdwpttchi2_user_chain(const char* path) {
    TChain* chain = new TChain("muons");
    // TChain* chain = new TChain("tree");
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
    unsigned int ntracks_cdclassify, ntracks_wpclassify;
    unsigned int nstoppings_cdclassify, nstoppings_wpclassify;

    chain->SetBranchAddress("run_id", &run_id);
    chain->SetBranchAddress("sec", &sec);
    chain->SetBranchAddress("nsec", &nsec);
    // chain->SetBranchAddress("totq_cd", &totq_cd);
    // chain->SetBranchAddress("totq_wp", &totq_wp);
    chain->SetBranchAddress("chi2", &chi2);
    chain->SetBranchAddress("iposx", &iposx);
    chain->SetBranchAddress("iposy", &iposy);
    chain->SetBranchAddress("iposz", &iposz);
    chain->SetBranchAddress("fposx", &fposx);
    chain->SetBranchAddress("fposy", &fposy);
    chain->SetBranchAddress("fposz", &fposz);
    // chain->SetBranchAddress("ntracks_cdclassify", &ntracks_cdclassify);
    // chain->SetBranchAddress("ntracks_wpclassify", &ntracks_wpclassify);
    // chain->SetBranchAddress("nstoppings_cdclassify", &nstoppings_cdclassify);
    // chain->SetBranchAddress("nstoppings_wpclassify", &nstoppings_wpclassify);

    long nentries = chain->GetEntries();
    std::cout << "Info: Found " << nentries << " entries in CdWpTtChi2 files\n";
    for (long k = 0l; k < nentries; ++k) {
        chain->GetEntry(k);
        tracks.insert(track{
            .run_id = run_id,
            .ts = TTimeStamp(sec, nsec),
            .totq_cd = 0.0, // totq_cd,
            .totq_wp = 0.0, // totq_wp,
            .quality = chi2,
            .ipos = TVector3(iposx, iposy, iposz),
            .fpos = TVector3(fposx, fposy, fposz),
            .is_single = true, // (ntracks_wpclassify == 1),
            .is_stopping = false // (nstoppings_wpclassify > 0)
        });
    }
    return tracks;
}

std::set<track> open_janus_user_chain(const char* path) {
    TChain* chain = new TChain("janus");
    chain->Add(path);
    std::set<track> tracks;

    int run_id;
    time_t sec;
    int nsec;
    double iposx, iposy, iposz;
    double fposx, fposy, fposz;

    chain->SetBranchAddress("run_id", &run_id);
    chain->SetBranchAddress("sec", &sec);
    chain->SetBranchAddress("nsec", &nsec);
    chain->SetBranchAddress("iposx", &iposx);
    chain->SetBranchAddress("iposy", &iposy);
    chain->SetBranchAddress("iposz", &iposz);
    chain->SetBranchAddress("fposx", &fposx);
    chain->SetBranchAddress("fposy", &fposy);
    chain->SetBranchAddress("fposz", &fposz);

    long nentries = chain->GetEntries();
    std::cout << "Info: Found " << nentries << " entries in Janus files\n";
    for (long k = 0l; k < nentries; ++k) {
        chain->GetEntry(k);
        tracks.insert(track{
            .run_id = run_id,
            .ts = TTimeStamp(sec, nsec),
            .totq_cd = 0.0,
            .totq_wp = 0.0,
            .quality = 0.0,
            .ipos = TVector3(iposx, iposy, iposz),
            .fpos = TVector3(fposx, fposy, fposz),
            .is_single = true,
            .is_stopping = false
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
        // if (!has_tt_info) continue;
        // if (!is_in_acrylic) continue; // SELECTION!
        if (ntracks_cdclassify != 1) continue; // SELECTION!
        if (ntracks_wpclassify != 1) continue; // SELECTION! 
        // if (stopping_cdclassify) continue; // SELECTION!
        // if (stopping_wpclassify) continue; // SELECTION!
        for (std::size_t i = 0ul; i < method->size(); ++i) {
            if ((*method)[i] == "CdWpTtChi2") continue;
            tracks[(*method)[i]].insert(track{
                .run_id = run_id,
                .ts = TTimeStamp(sec, nsec),
                .totq_cd = totq_cd,
                .totq_wp = totq_wp,
                .quality = (*quality)[i],
                .ipos = TVector3((*iposx)[i], (*iposy)[i], (*iposz)[i]),
                .fpos = TVector3((*fposx)[i], (*fposy)[i], (*fposz)[i]),
                .is_single = (ntracks_wpclassify == 1),
                .is_stopping = stopping_wpclassify
            });
        }
    }
    return tracks;
}

struct MuonPerformance {
    double angle;
    double distance;
    double clippingness;
    double clippingness_trk;
    int run_id;
    time_t sec;
    int nsec;
    double quality;
    double tt_quality;
    double zenith;
    double azimuth;
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
        performances[method] = {};

        std::cout << "\n--- Correlating " << method << " with Tt ---" << std::endl;

        for (const track& trk : track_set) {
            TTimeStamp lower_bound_ts(trk.ts.GetSec(), trk.ts.GetNanoSec() - 1000);
            TTimeStamp upper_bound_ts(trk.ts.GetSec(), trk.ts.GetNanoSec() + 1000);
            
            std::set<track>::const_iterator it_tt = tt_tracks.lower_bound({0, lower_bound_ts, 0, 0, {}, {}});

            while (it_tt != tt_tracks.end() && lower_bound_ts <= it_tt->ts && it_tt->ts <= upper_bound_ts) {
                performances[method].push_back(MuonPerformance{
                    .angle = compute_angle_between_track(trk, *it_tt),
                    .distance = compute_distance_between_track(trk, *it_tt),
                    .clippingness = compute_clippingness(*it_tt),
                    .clippingness_trk = compute_clippingness(trk),
                    .run_id = trk.run_id,
                    .sec = trk.ts.GetSec(),
                    .nsec = trk.ts.GetNanoSec(),
                    .quality = trk.quality,
                    .tt_quality = it_tt->quality,
                    .zenith = (trk.fpos - trk.ipos).Unit().Theta() * 180.0 / M_PI,
                    .azimuth = (trk.fpos - trk.ipos).Unit().Phi() * 180.0 / M_PI
                });
                ++it_tt;
            }
        }
    }
    return performances;
}

std::map<std::string, std::vector<MuonPerformance>> compute_global_correlations(std::map<std::string, std::set<track>>& tracks, const std::string& refname) {
    if (tracks.find(refname) == tracks.end()) {
        std::cerr << "Error: " << refname << " tracks not found in map.\n";
        return {};
    }
    const std::set<track>& ref_tracks = tracks[refname];

    std::map<std::string, std::vector<MuonPerformance>> performances;

    for (const track& ref_muon : ref_tracks) {
        std::map<std::string, track> coincident_map;
        bool all_found = true;

        for (const auto& [method, track_set] : tracks) {
            if (method == refname) continue;
            bool found_in_method = false;
            TTimeStamp lower_bound_ts(ref_muon.ts.GetSec(), ref_muon.ts.GetNanoSec() - 1000);
            TTimeStamp upper_bound_ts(ref_muon.ts.GetSec(), ref_muon.ts.GetNanoSec() + 1000);
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
                    .angle = compute_angle_between_track(ref_muon, muon),
                    .distance = compute_distance_between_track(ref_muon, muon),
                    .clippingness = compute_clippingness(ref_muon),
                    .clippingness_trk = compute_clippingness(muon),
                    .run_id = muon.run_id,
                    .sec = muon.ts.GetSec(),
                    .nsec = muon.ts.GetNanoSec(),
                    .quality = muon.quality,
                    .tt_quality = ref_muon.quality,
                    .zenith = (muon.fpos - muon.ipos).Unit().Theta() * 180.0 / M_PI,
                    .azimuth = (muon.fpos - muon.ipos).Unit().Phi() * 180.0 / M_PI
                });
            }
        }
    }
    return performances;
}

std::map<std::string, std::vector<MuonPerformance>> compute_global_correlations_no_tt(std::map<std::string, std::set<track>>& tracks, const std::string& refname) {
    if (tracks.find(refname) == tracks.end()) {
        std::cerr << "Error: " << refname << " tracks not found in map.\n";
        return {};
    }
    const std::set<track>& ref_tracks = tracks[refname];

    std::map<std::string, std::vector<MuonPerformance>> performances;

    for (const track& ref_muon : ref_tracks) {
        std::map<std::string, track> coincident_map;
        bool all_found = true;

        for (const auto& [method, track_set] : tracks) {
            if (method == "Tt") continue;
            bool found_in_method = false;
            TTimeStamp lower_bound_ts(ref_muon.ts.GetSec(), ref_muon.ts.GetNanoSec() - 1000);
            TTimeStamp upper_bound_ts(ref_muon.ts.GetSec(), ref_muon.ts.GetNanoSec() + 1000);
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
                    .angle = compute_angle_between_track(ref_muon, muon),
                    .distance = compute_distance_between_track(ref_muon, muon),
                    .clippingness = compute_clippingness(ref_muon),
                    .clippingness_trk = compute_clippingness(muon),
                    .run_id = muon.run_id,
                    .sec = muon.ts.GetSec(),
                    .nsec = muon.ts.GetNanoSec(),
                    .quality = muon.quality,
                    .tt_quality = ref_muon.quality,
                    .zenith = (muon.fpos - muon.ipos).Unit().Theta() * 180.0 / M_PI,
                    .azimuth = (muon.fpos - muon.ipos).Unit().Phi() * 180.0 / M_PI
                });
            }
        }
    }
    return performances;
}

int fast_muon_reconstruction_comparison(const char* path_joint, const char* path_cdwpttchi2, const char* path_amber, const char* path_edwin, const char* path_janus) {
    std::map<std::string, std::set<track>> joint_tracks = open_joint_reco_user_chain(path_joint);
    std::set<track> cdclassify_tracks = joint_tracks["CdClassify"];
    std::set<track> wpclassify_tracks = joint_tracks["WpBasic"];
    std::set<track> tt_tracks = joint_tracks["Tt"];
    std::set<track> cdwpttchi2_tracks = open_cdwpttchi2_user_chain(path_cdwpttchi2);
    std::set<track> amber_v5_5_tracks = open_amber_v5_5_user_chain(path_amber);
    std::set<track> edwin_tracks = open_edwin_user_chain(path_edwin);
    std::set<track> janus_tracks = open_janus_user_chain(path_janus);

    std::map<std::string, std::set<track>> tracks;
    tracks["Tt"] = tt_tracks;
    tracks["CdClassify"] = cdclassify_tracks;
    tracks["WpBasic"] = wpclassify_tracks;
    if (!cdwpttchi2_tracks.empty()) {
        tracks["CdWpTtChi2"] = cdwpttchi2_tracks;
    } 
    if (!amber_v5_5_tracks.empty()) {
        tracks["Amber_v5.5"] = amber_v5_5_tracks;
    }
    if (!edwin_tracks.empty()) {
        tracks["Edwin"] = edwin_tracks;
    }
    if (!janus_tracks.empty()) {
        tracks["Janus"] = janus_tracks;
    }

    for (const auto& [method, track_set] : tracks) {
        std::cout << method << " size: " << track_set.size() << '\n';
    }

    // std::map<std::string, std::vector<MuonPerformance>> performances = compute_correlations(tracks);
    std::map<std::string, std::vector<MuonPerformance>> performances = compute_global_correlations(tracks, "Tt");
    std::map<std::string, std::vector<MuonPerformance>> performances_no_tt = compute_global_correlations_no_tt(tracks, "CdWpTtChi2");
    
    std::map<std::string, std::vector<double>> angles;
    std::map<std::string, std::vector<double>> distances;
    for (const auto& [method, perf] : performances) {
        std::cout << "Extracting angles and distances for " << method << " (size = " << perf.size() << ")\n";
        angles[method] = extract_angles_from_performances(perf);
        distances[method] = extract_distances_from_performances(perf);
    }

    std::map<std::string, TH1D*> method_angle_map;
    std::map<std::string, TH1D*> method_distance_map;

    // double xmin_angle = 0.0, xmax_angle = 180.0;
    // double xmin_distance = 0.0, xmax_distance = 40.0;
    // int nbins_angle = 500, nbins_distance = 500;
    double xmin_angle = 0.0, xmax_angle = 5.0;
    double xmin_distance = 0.0, xmax_distance = 2.0;
    int nbins_angle = 50, nbins_distance = 50;

    method_angle_map["CdClassify"] = new TH1D("h_angle_cdclassify", "Angle between tracks direction (CdClassify);#alpha (deg);Entries;", nbins_angle, xmin_angle, xmax_angle);
    method_distance_map["CdClassify"] = new TH1D("h_distance_cdclassify", "Distance between tracks middle point (CdClassify);d_{mid} (m);Entries;", nbins_distance, xmin_distance, xmax_distance);
    method_angle_map["WpBasic"] = new TH1D("h_angle_wpclassify", "Angle between tracks direction (WpClassify);#alpha (deg);Entries;", nbins_angle, xmin_angle, xmax_angle);
    method_distance_map["WpBasic"] = new TH1D("h_distance_wpclassify", "Distance between tracks middle point (WpClassify);d_{mid} (m);Entries;", nbins_distance, xmin_distance, xmax_distance);
    if (!cdwpttchi2_tracks.empty()) {
        method_angle_map["CdWpTtChi2"] = new TH1D("h_angle_cdwpttchi2", "Angle between tracks direction (CdWpTtChi2);#alpha (deg);Entries;", nbins_angle, xmin_angle, xmax_angle);
        method_distance_map["CdWpTtChi2"] = new TH1D("h_distance_cdwpttchi2", "Distance between tracks middle point (CdWpTtChi2);d_{mid} (m);Entries;", nbins_distance, xmin_distance, xmax_distance);
    }
    if (!amber_v5_5_tracks.empty()) {
        method_angle_map["Amber_v5.5"] = new TH1D("h_angle_amber", "Angle between tracks direction (Amber);#alpha (deg);Entries;", nbins_angle, xmin_angle, xmax_angle);
        method_distance_map["Amber_v5.5"] = new TH1D("h_distance_amber", "Distance between tracks middle point (Amber);d_{mid} (m);Entries;", nbins_distance, xmin_distance, xmax_distance);
    }
    if (!edwin_tracks.empty()) {
        method_angle_map["Edwin"] = new TH1D("h_angle_edwin", "Angle between tracks direction (Edwin);#alpha (deg);Entries;", nbins_angle, xmin_angle, xmax_angle);
        method_distance_map["Edwin"] = new TH1D("h_distance_edwin", "Distance between tracks middle point (Edwin);d_{mid} (m);Entries;", nbins_distance, xmin_distance, xmax_distance);
    }
    if (!janus_tracks.empty()) {
        method_angle_map["Janus"] = new TH1D("h_angle_janus", "Angle between tracks direction (Janus);#alpha (deg);Entries;", nbins_angle, xmin_angle, xmax_angle);
        method_distance_map["Janus"] = new TH1D("h_distance_janus", "Distance between tracks middle point (Janus);d_{mid} (m);Entries;", nbins_distance, xmin_distance, xmax_distance);
    }


    TFile* fout = TFile::Open("/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/performance/output.root", "RECREATE");
    if (!fout) {
        std::cerr << "Cannot open output file output.root\n";
        return 1;
    }
    fout->cd();

    int run_id;
    time_t sec;
    int nsec;
    double angle;
    double chi2;
    double dist_center;
    double dist_mid_point;
    double zenith;
    double azimuth;

    std::map<std::string, TTree*> trees;
    for (const auto& [method, perf] : performances) {
        trees[method] = new TTree(method.c_str(), method.c_str());
        trees[method]->Branch("run_id", &run_id);
        trees[method]->Branch("sec", &sec);
        trees[method]->Branch("nsec", &nsec);
        trees[method]->Branch("angle", &angle);
        trees[method]->Branch("chi2", &chi2);
        trees[method]->Branch("dist_center", &dist_center);
        trees[method]->Branch("dist_mid_point", &dist_mid_point);
        trees[method]->Branch("zenith", &zenith);
        trees[method]->Branch("azimuth", &azimuth);

        for (const MuonPerformance& mp : perf) {
            method_angle_map[method]->Fill(mp.angle);
            method_distance_map[method]->Fill(mp.distance);
            run_id = mp.run_id;
            sec = mp.sec;
            nsec = mp.nsec;
            angle = mp.angle;
            chi2 = mp.quality;
            dist_center = mp.clippingness;
            dist_mid_point = mp.distance;
            zenith = mp.zenith;
            azimuth = mp.azimuth;
            trees[method]->Fill();
        }
    }

    fout->cd();
    for (const auto& [method, t] : trees) {
        t->Write();
    }
    fout->Close();

    TFile* fout_no_tt = TFile::Open("/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/performance/output_no_tt.root", "RECREATE");
    if (!fout_no_tt) {
        std::cerr << "Cannot open output file output_no_tt.root\n";
        return 1;
    }
    fout_no_tt->cd();

    std::map<std::string, TTree*> trees_no_tt;
    for (const auto& [method, perf] : performances_no_tt) {
        trees_no_tt[method] = new TTree((method + "_no_tt").c_str(), (method + "_no_tt").c_str());
        trees_no_tt[method]->Branch("run_id", &run_id);
        trees_no_tt[method]->Branch("sec", &sec);
        trees_no_tt[method]->Branch("nsec", &nsec);
        trees_no_tt[method]->Branch("angle", &angle);
        trees_no_tt[method]->Branch("chi2", &chi2);
        trees_no_tt[method]->Branch("dist_center", &dist_center);
        trees_no_tt[method]->Branch("dist_mid_point", &dist_mid_point);
        trees_no_tt[method]->Branch("zenith", &zenith);
        trees_no_tt[method]->Branch("azimuth", &azimuth);

        std::cout << "Size performance no TT (" << method << "): " << perf.size() << '\n';

        for (const MuonPerformance& mp : perf) {
            run_id = mp.run_id;
            sec = mp.sec;
            nsec = mp.nsec;
            angle = mp.angle;
            chi2 = mp.quality;
            dist_center = mp.clippingness;
            dist_mid_point = mp.distance;
            zenith = mp.zenith;
            azimuth = mp.azimuth;
            trees_no_tt[method]->Fill();
        }
    }

    fout_no_tt->cd();
    for (const auto& [method, t] : trees_no_tt) {
        t->Write();
    }
    fout_no_tt->Close();

    return 0;
}