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

double get_quantile(std::vector<double>::const_iterator first, std::vector<double>::const_iterator last, double quantile) {
    if (first == last) return 0.0;
    std::vector<double> v(first, last);
    std::vector<double>::iterator position = v.begin() + quantile * v.size();
    std::nth_element(v.begin(), position, v.end());
    return *position;
}

void plot_metrics(const std::vector<TH1D*>& hists, const std::vector<double>& quantiles) {
    if (hists.empty()) return;
    std::vector<Color_t> colors = {kBlack, kGreen + 2, kViolet, kBlue, kRed};
    TCanvas* c = new TCanvas(Form("c_%s", hists[0]->GetName()), "Metric", 1000, 1000);
    c->cd();

    double max = 0.0;
    for (std::size_t i = 0ul; i < hists.size(); ++i) {
        TH1D* h = hists[i];
        if (h->GetMaximum() > max) {
            max = h->GetMaximum();
        }
        std::cout << "68.2% " << h->GetName() << ": " << quantiles[i] << '\n';
    }
    hists[0]->SetMaximum(max * 1.1);

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
        TLine* line = new TLine(quantiles[i], 0.0, quantiles[i], h->GetMaximum());
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->SetLineColor(colors[i]);
        line->Draw();
    };

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
    std::cout << "Info: Found " << nentries << " entries in Amber_v5.5 files\n";
    for (long k = 0l; k < nentries; ++k) {
        chain->GetEntry(k);
        tracks.insert(track{
            .run_id = runID,
            .ts = TTimeStamp(fSec, fNanoSec),
            .totq_cd = 0.0,
            .totq_wp = charge,
            .ipos = TVector3(xin, yin, zin),
            .fpos = TVector3(xout, yout, zout)
        });
    }
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
    std::cout << "Info: Found " << nentries << " entries in EDWIN files\n";
    for (long k = 0l; k < nentries; ++k) {
        chain->GetEntry(k);
        tracks.insert(track{
            .run_id = run_number,
            .ts = TTimeStamp(cd_time_s, cd_time_ns),
            .totq_cd = cd_totalPE,
            .totq_wp = wp_totalPE,
            .ipos = TVector3(enterX, enterY, enterZ),
            .fpos = TVector3(exitX, exitY, exitZ)
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
        int ntracks_cdclassify = 0;
        int ntracks_wpclassify = 0;
        for (std::size_t i = 0ul; i < method->size(); ++i) {
            if ((*method)[i] == "Tt") {
                has_tt_info = true;
            }
            if ((*method)[i] == "CdClassify") {
                ++ntracks_cdclassify;
            }
            if ((*method)[i] == "WpBasic") {
                ++ntracks_wpclassify;
            }
        }
        if (!has_tt_info) continue;
        if (ntracks_cdclassify != 1) continue;
        if (ntracks_wpclassify != 1) continue;
        for (std::size_t i = 0ul; i < method->size(); ++i) {
            tracks[(*method)[i]].insert(track{
                .run_id = run_id,
                .ts = TTimeStamp(sec, nsec),
                .totq_cd = totq_cd,
                .totq_wp = totq_wp,
                .ipos = TVector3((*iposx)[i], (*iposy)[i], (*iposz)[i]),
                .fpos = TVector3((*fposx)[i], (*fposy)[i], (*fposz)[i])
            });
        }
    }
    return tracks;
}

std::pair<std::map<std::string, std::vector<double>>, std::map<std::string, std::vector<double>>> compute_correlations(std::map<std::string, std::set<track>>& tracks) {
    if (tracks.find("Tt") == tracks.end()) {
        std::cerr << "Error: Tt tracks not found in map.\n";
        return {};
    }
    const std::set<track>& tt_tracks = tracks["Tt"];

    std::map<std::string, std::vector<double>> angles;
    std::map<std::string, std::vector<double>> distances;

    for (const auto& [method, track_set] : tracks) {
        if (method == "Tt") continue;

        std::cout << "\n--- Correlating " << method << " with Tt ---" << std::endl;

        for (const track& trk : track_set) {
            TTimeStamp lower_bound_ts(trk.ts.GetSec(), trk.ts.GetNanoSec() - 1000);
            TTimeStamp upper_bound_ts(trk.ts.GetSec(), trk.ts.GetNanoSec() + 1000);
            
            std::set<track>::const_iterator it_tt = tt_tracks.lower_bound({0, lower_bound_ts, 0, 0, {}, {}});

            while (it_tt != track_set.end() && lower_bound_ts <= it_tt->ts && it_tt->ts <= upper_bound_ts) {
                angles[method].push_back(compute_angle_between_track(trk, *it_tt));
                distances[method].push_back(compute_distance_between_track(trk, *it_tt));
                ++it_tt;
            }
        }
    }
    return std::make_pair(angles, distances);
}

std::pair<std::map<std::string, std::vector<double>>, std::map<std::string, std::vector<double>>> compute_global_correlations(std::map<std::string, std::set<track>>& tracks) {
    if (tracks.find("Tt") == tracks.end()) {
        std::cerr << "Error: Tt tracks not found in map.\n";
        return {};
    }
    const std::set<track>& tt_tracks = tracks["Tt"];

    std::map<std::string, std::vector<double>> angles;
    std::map<std::string, std::vector<double>> distances;

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
                angles[method].push_back(compute_angle_between_track(muon, tt_muon));
                distances[method].push_back(compute_distance_between_track(muon, tt_muon));
            }
        }
    }
    return std::make_pair(angles, distances);
}

int fast_muon_reconstruction_comparison(const char* path_joint, const char* path_amber, const char* path_edwin) {
    std::map<std::string, std::set<track>> tracks = open_joint_reco_user_chain(path_joint);
    tracks["Amber_v5.5"] = open_amber_v5_5_user_chain(path_amber);
    tracks["Edwin"] = open_edwin_user_chain(path_edwin);

    std::pair<std::map<std::string, std::vector<double>>, std::map<std::string, std::vector<double>>> correlations = compute_correlations(tracks);
    const std::map<std::string, std::vector<double>>& angles = correlations.first;
    const std::map<std::string, std::vector<double>>& distances = correlations.second;

    std::map<std::string, TH1D*> method_angle_map;
    std::map<std::string, TH1D*> method_distance_map;

    method_angle_map["CdWpTtChi2"] = new TH1D("h_angle_cdwpttchi2", "Angle between tracks direction (CdWpTtChi2);#alpha (deg);Entries;", 25, 0.0, 5.0);
    method_distance_map["CdWpTtChi2"] = new TH1D("h_distance_cdwpttchi2", "Distance between tracks middle point (CdWpTtChi2);d_{mid} (m);Entries;", 25, 0.0, 2.0);
    method_angle_map["CdClassify"] = new TH1D("h_angle_cdclassify", "Angle between tracks direction (CdClassify);#alpha (deg);Entries;", 25, 0.0, 5.0);
    method_distance_map["CdClassify"] = new TH1D("h_distance_cdclassify", "Distance between tracks middle point (CdClassify);d_{mid} (m);Entries;", 25, 0.0, 2.0);
    method_angle_map["WpBasic"] = new TH1D("h_angle_wpclassify", "Angle between tracks direction (WpClassify);#alpha (deg);Entries;", 25, 0.0, 5.0);
    method_distance_map["WpBasic"] = new TH1D("h_distance_wpclassify", "Distance between tracks middle point (WpClassify);d_{mid} (m);Entries;", 25, 0.0, 2.0);
    method_angle_map["Amber_v5.5"] = new TH1D("h_angle_amber", "Angle between tracks direction (Amber);#alpha (deg);Entries;", 25, 0.0, 5.0);
    method_distance_map["Amber_v5.5"] = new TH1D("h_distance_amber", "Distance between tracks middle point (Amber);d_{mid} (m);Entries;", 25, 0.0, 2.0);
    method_angle_map["Edwin"] = new TH1D("h_angle_edwin", "Angle between tracks direction (Edwin);#alpha (deg);Entries;", 25, 0.0, 5.0);
    method_distance_map["Edwin"] = new TH1D("h_distance_edwin", "Distance between tracks middle point (Edwin);d_{mid} (m);Entries;", 25, 0.0, 2.0);

    for (const auto& [method, ang] : angles) {
        for (std::vector<double>::const_iterator it = ang.begin(); it != ang.end(); ++it) {
            method_angle_map[method]->Fill(*it);
        }
    }
    for (const auto& [method, dist] : distances) {
        for (std::vector<double>::const_iterator it = dist.begin(); it != dist.end(); ++it) {
            method_distance_map[method]->Fill(*it);
        }
    }

    // std::cout << "68.2% angle: " << get_quantile(angles.begin(), angles.end(), 0.682) << '\n';
    // std::cout << "95.4% angle: " << get_quantile(angles.begin(), angles.end(), 0.954) << '\n';
    // std::cout << "68.2% distance: " << get_quantile(distances.begin(), distances.end(), 0.682) << '\n';
    // std::cout << "95.4% distance: " << get_quantile(distances.begin(), distances.end(), 0.954) << '\n';

    std::vector<double> angle_quantiles = {
        get_quantile(method_angle_map["CdWpTtChi2"].begin(), method_angle_map["CdWpTtChi2"].end(), 0.682),
        get_quantile(method_angle_map["CdClassify"].begin(), method_angle_map["CdClassify"].end(), 0.682),
        get_quantile(method_angle_map["WpBasic"].begin(), method_angle_map["WpBasic"].end(), 0.682),
        get_quantile(method_angle_map["Amber_v5.5"].begin(), method_angle_map["Amber_v5.5"].end(), 0.682),
        get_quantile(method_angle_map["Edwin"].begin(), method_angle_map["Edwin"].end(), 0.682)
    };
    std::vector<double> distance_quantiles = {
        get_quantile(method_distance_map["CdWpTtChi2"].begin(), method_distance_map["CdWpTtChi2"].end(), 0.682),
        get_quantile(method_distance_map["CdClassify"].begin(), method_distance_map["CdClassify"].end(), 0.682),
        get_quantile(method_distance_map["WpBasic"].begin(), method_distance_map["WpBasic"].end(), 0.682),
        get_quantile(method_distance_map["Amber_v5.5"].begin(), method_distance_map["Amber_v5.5"].end(), 0.682),
        get_quantile(method_distance_map["Edwin"].begin(), method_distance_map["Edwin"].end(), 0.682)
    };


    plot_metrics({method_angle_map["CdWpTtChi2"], method_angle_map["CdClassify"], method_angle_map["WpBasic"], method_angle_map["Amber_v5.5"], method_angle_map["Edwin"]}, angle_quantiles);
    plot_metrics({method_distance_map["CdWpTtChi2"], method_distance_map["CdClassify"], method_distance_map["WpBasic"], method_distance_map["Amber_v5.5"], method_distance_map["Edwin"]}, distance_quantiles);

    return 0;
}