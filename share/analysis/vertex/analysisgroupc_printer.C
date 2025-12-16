#include <algorithm>
#include <chrono>
#include <iostream>
#include <set>
#include <string>

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include "analysis.hpp"
#include "utils/event.hpp"
#include "utils/timestamp.hpp"
#include "utils/vec3.hpp"

struct IBD {

    vertex prompt;
    vertex delayed;

};

struct Cosmo {

    vertex prompt;
    vertex delayed;
    double dlat_mu2p;
    double dlat_mu2d;
    double dt_mu2p;
    double dt_mu2d;

};

inline bool operator<(const IBD& lhs, const IBD& rhs) {
    return lhs.prompt.ts < rhs.prompt.ts;
}

inline bool operator<(const Cosmo& lhs, const Cosmo& rhs) {
    return lhs.prompt.ts < rhs.prompt.ts;
}

std::vector<IBD> get_all_ibd(const std::string& filename, analysis_base* analysis) {
    if (!analysis->retrieve(filename)) return {};
    std::set<IBD> ibds_ordered;
    std::cout << "=== Analysis: " << analysis->name << " (Total Entries: " << analysis->size() << ") ===\n";
    for (long k = 0; k < analysis->size(); ++k) {
        analysis->entry(k);
        if (k % 1000 == 0) {
            std::cout << "\rProcessing Entry " << k << " / " << analysis->size();
        }
        // analysis->print();
        if (!analysis->selection()) continue;
        ibds_ordered.insert({analysis->prompt, analysis->delayed});
    }

    std::vector<IBD> ibds;
    ibds.reserve(ibds_ordered.size());
    for (std::set<IBD>::const_iterator it = ibds_ordered.begin(); it != ibds_ordered.end(); ++it) {
        ibds.push_back(*it);
    }
    return ibds;
}

std::vector<Cosmo> get_all_cosmo(const std::string& filename, cosmo_shape_analysis* analysis) {
    if (!analysis->retrieve(filename)) return {};
    std::set<Cosmo> cosmos_ordered;
    std::cout << "=== Analysis: Cosmo (Total Entries: " << analysis->size() << ") ===\n";
    for (long k = 0; k < analysis->size(); ++k) {
        analysis->entry(k);
        if (k % 1000 == 0) {
            std::cout << "\rProcessing Entry " << k << " / " << analysis->size();
        }
        // analysis->print();
        if (!analysis->selection()) continue;
        Cosmo cosmo;
        cosmo.prompt = analysis->prompt;
        cosmo.delayed = analysis->delayed;
        cosmo.dlat_mu2p = analysis->dlat_mu2p;
        cosmo.dlat_mu2d = analysis->dlat_mu2d;
        cosmo.dt_mu2p = analysis->dt_mu2p;
        cosmo.dt_mu2d = analysis->dt_mu2d;
        cosmos_ordered.insert(cosmo);
    }

    std::vector<Cosmo> cosmos;
    cosmos.reserve(cosmos_ordered.size());
    for (std::set<Cosmo>::const_iterator it = cosmos_ordered.begin(); it != cosmos_ordered.end(); ++it) {
        cosmos.push_back(*it);
    }
    return cosmos;
}

void fit_and_plot_cosmo_rate_with_neutron(TH1D* h) {
    TCanvas* c = new TCanvas(Form("%s_canvas", h->GetName()), Form("%s Canvas", h->GetName()), 1000, 1000);
    c->cd();

    double constant_term = 0.0;
    for (int bin = h->GetXaxis()->FindBin(0.8); bin <= h->GetXaxis()->GetNbins(); ++bin) {
        constant_term += h->GetBinContent(bin);
    }
    constant_term /= (h->GetXaxis()->GetNbins() - h->GetXaxis()->FindBin(0.8) + 1);
    double exponential_term = h->GetMaximum() - constant_term;

    TF1* f = new TF1(Form("f_%s", h->GetName()), "[0] + [1] * exp(-x / [2])", 0.02, 1.2);
    f->SetParameter(0, constant_term);
    f->SetParameter(1, exponential_term);
    f->SetParameter(2, 180.0e-3);

    f->SetLineColor(kRed);
    f->SetLineWidth(3);

    h->Fit(f, "R");
    h->SetLineWidth(3);
    h->GetXaxis()->SetTitle("#Delta t_{#mu2p} (s)");
    h->GetXaxis()->CenterTitle(kTRUE);
    h->GetYaxis()->SetTitle("Entries");
    h->GetYaxis()->CenterTitle(kTRUE);
    h->GetYaxis()->SetTitleOffset(1.5);
    h->Draw();
    f->Draw("SAME");
    
    c->SetTickx();
    c->SetTicky();

    c->Update();

    double time_window = 1.2;
    double binning = time_window / 120.0;
    std::cout << "Fit Results for " << h->GetName() << ":\n";
    std::cout << "nIBD = " << f->GetParameter(0) * time_window / binning << " +/- " << f->GetParError(0) * time_window / binning << '\n';
    std::cout << "nLiHe = " << f->GetParameter(1) * f->GetParameter(2) * (1 - std::exp(-time_window / f->GetParameter(2))) / binning << '\n';    
}

struct PhysicalMuon {

    timestamp ts;
    std::vector<std::size_t> indices;

};

struct MuonDataAssociation {

    double dt;
    int neutron_count;
    std::vector<double> dlat_cdwp;
    std::vector<double> dlat_tt;

};

void analyze_cosmo_rate_with_neutron(const std::string& filename, cosmo_rate_analysis* analysis) {
    if (!analysis->retrieve(filename)) return;
    std::map<IBD, std::vector<MuonDataAssociation>> ibds_to_mu;
    std::cout << "=== Analysis: CosmoRateWithNeutron (Total Entries: " << analysis->size() << ") ===\n";
    for (long k = 0; k < analysis->size(); ++k) {
        analysis->entry(k);
        if (k % 1000 == 0) {
            std::cout << "\rProcessing Entry " << k << " / " << analysis->size();
        }
        // analysis->print();
        if (!analysis->selection()) continue;
        IBD ibd;
        ibd.prompt = analysis->prompt;
        ibd.delayed = analysis->delayed;

        std::vector<PhysicalMuon> physical_muon;
        std::vector<MuonDataAssociation> muon_data;
        for (std::size_t k = 0ul; k < analysis->method_mu->size(); ++k) {
            timestamp ts_mu{analysis->sec_mu->operator[](k), analysis->nsec_mu->operator[](k)};
            std::vector<PhysicalMuon>::iterator it = std::find_if(
                physical_muon.begin(), physical_muon.end(),
                [ts_mu](const PhysicalMuon& mu) {
                    timestamp diff = ts_mu - mu.ts;
                    return timestamp{0, -1000} < diff && diff < timestamp{0, 1000};
                }
            );
            if (it == physical_muon.end()) {
                physical_muon.push_back({ts_mu, {k}});
            }
            else {
                it->indices.push_back(k);
            }
        }

        for (const PhysicalMuon& mu : physical_muon) {
            timestamp diff = ibd.prompt.ts - mu.ts;
            if (diff < timestamp{0, 5000000} || timestamp{0, 1200000000} < diff) continue;

            MuonDataAssociation assoc;
            assoc.dt = timestamp_to_double(diff);

            int neutron_count = 0;
            for (std::size_t k = 0ul; k < analysis->sec_n->size(); ++k) {
                double e_n = analysis->e_n->operator[](k);
                if (e_n < 2.0 || 2.5 < e_n) continue;
                timestamp ts_n{analysis->sec_n->operator[](k), analysis->nsec_n->operator[](k)};
                if (ts_n < mu.ts + timestamp{0, 20000} || mu.ts + timestamp{0, 2000000} < ts_n) continue;
                ++neutron_count;
            }
            assoc.neutron_count = neutron_count;

            for (std::size_t idx : mu.indices) {
                const std::string& method = analysis->method_mu->operator[](idx);
                if (method != "CdWpTtChi2" && method != "Tt") continue;
                vec3 pos_mu{analysis->posx_mu->operator[](idx), analysis->posy_mu->operator[](idx), analysis->posz_mu->operator[](idx)};
                vec3 dir_mu{analysis->dirx_mu->operator[](idx), analysis->diry_mu->operator[](idx), analysis->dirz_mu->operator[](idx)};
                double d = mag(cross(dir_mu, ibd.prompt.pos - pos_mu));

                if (method == "CdWpTtChi2") assoc.dlat_cdwp.push_back(d);
                else assoc.dlat_tt.push_back(d);
            }
            muon_data.push_back(std::move(assoc));
        }
        ibds_to_mu[ibd] = muon_data;
    }

    TH1D* h_cosmo_rate_with_neutron = new TH1D("h_cosmo_rate_with_neutron", "Cosmo Rate With Neutron", 120, 0.0, 1.2);
    TH1D* h_cosmo_rate_with_no_neutron = new TH1D("h_cosmo_rate_with_no_neutron", "Cosmo Rate With No Neutron", 120, 0.0, 1.2);
    TH1D* h_cosmo_rate_with_at_least_1_neutron = new TH1D("h_cosmo_rate_with_at_least_1_neutron", "Cosmo Rate With At Least 1 Neutron", 120, 0.0, 1.2);
    TH1D* h_cosmo_rate_with_at_least_2_neutron = new TH1D("h_cosmo_rate_with_at_least_2_neutron", "Cosmo Rate With At Least 2 Neutron", 120, 0.0, 1.2);
    TH1D* h_cosmo_rate_with_at_least_3_neutron = new TH1D("h_cosmo_rate_with_at_least_3_neutron", "Cosmo Rate With At Least 3 Neutron", 120, 0.0, 1.2);
    for (const std::pair<IBD, std::vector<MuonDataAssociation>>& val : ibds_to_mu) {
        const std::vector<MuonDataAssociation>& muon_data = val.second;
        for (const MuonDataAssociation& assoc : muon_data) {
            h_cosmo_rate_with_neutron->Fill(assoc.dt);
            if (assoc.neutron_count == 0) {
                h_cosmo_rate_with_no_neutron->Fill(assoc.dt);
            }
            if (assoc.neutron_count >= 1) {
                h_cosmo_rate_with_at_least_1_neutron->Fill(assoc.dt);
            }
            if (assoc.neutron_count >= 2) {
                h_cosmo_rate_with_at_least_2_neutron->Fill(assoc.dt);
            }
            if (assoc.neutron_count >= 3) {
                h_cosmo_rate_with_at_least_3_neutron->Fill(assoc.dt);
            }
        }
    }

    fit_and_plot_cosmo_rate_with_neutron(h_cosmo_rate_with_neutron);
    fit_and_plot_cosmo_rate_with_neutron(h_cosmo_rate_with_no_neutron);
    fit_and_plot_cosmo_rate_with_neutron(h_cosmo_rate_with_at_least_1_neutron);
    fit_and_plot_cosmo_rate_with_neutron(h_cosmo_rate_with_at_least_2_neutron);
    fit_and_plot_cosmo_rate_with_neutron(h_cosmo_rate_with_at_least_3_neutron);

    TH2D* h_d_mu2p_cdwp_vs_dt_mu2p = new TH2D("h_d_mu2p_cdwp_vs_dt_mu2p", "Cosmo time vs distance", 120, 0.0, 1.5, 100, 0.0, 40000.0);
    TH2D* h_d_mu2p_tt_vs_dt_mu2p = new TH2D("h_d_mu2p_tt_vs_dt_mu2p", "Cosmo time vs distance", 120, 0.0, 1.5, 100, 0.0, 40000.0);
    for (const std::pair<IBD, std::vector<MuonDataAssociation>>& val : ibds_to_mu) {
        const std::vector<MuonDataAssociation>& muon_data = val.second;
        for (const MuonDataAssociation& assoc : muon_data) {
            std::vector<double>::const_iterator it_dlat_cdwp = std::min_element(assoc.dlat_cdwp.begin(), assoc.dlat_cdwp.end());
            std::vector<double>::const_iterator it_dlat_tt = std::min_element(assoc.dlat_tt.begin(), assoc.dlat_tt.end());
            if (it_dlat_cdwp != assoc.dlat_cdwp.end()) {
                h_d_mu2p_cdwp_vs_dt_mu2p->Fill(assoc.dt, *it_dlat_cdwp);
            }
            if (it_dlat_tt != assoc.dlat_tt.end()) {
                h_d_mu2p_tt_vs_dt_mu2p->Fill(assoc.dt, *it_dlat_tt);
            }
        }
    }

    TCanvas* c_d_mu2p_cdwp_vs_dt_mu2p = new TCanvas("c_d_mu2p_cdwp_vs_dt_mu2p", "d_mu2p_cdwp vs dt_mu2p", 1000, 1000);
    c_d_mu2p_cdwp_vs_dt_mu2p->cd();
    h_d_mu2p_cdwp_vs_dt_mu2p->GetXaxis()->SetTitle("#Delta t_{#mu2p} (s)");
    h_d_mu2p_cdwp_vs_dt_mu2p->GetXaxis()->CenterTitle(kTRUE);
    h_d_mu2p_cdwp_vs_dt_mu2p->GetYaxis()->SetTitle("d_{#mu2p} (cm)");
    h_d_mu2p_cdwp_vs_dt_mu2p->GetYaxis()->CenterTitle(kTRUE);
    h_d_mu2p_cdwp_vs_dt_mu2p->GetYaxis()->SetTitleOffset(1.5);
    h_d_mu2p_cdwp_vs_dt_mu2p->Draw("COLZ");
    c_d_mu2p_cdwp_vs_dt_mu2p->SetTickx();
    c_d_mu2p_cdwp_vs_dt_mu2p->SetTicky();
    c_d_mu2p_cdwp_vs_dt_mu2p->Update();

    TCanvas* c_d_mu2p_tt_vs_dt_mu2p = new TCanvas("c_d_mu2p_tt_vs_dt_mu2p", "d_mu2p_tt vs dt_mu2p", 1000, 1000);
    c_d_mu2p_tt_vs_dt_mu2p->cd();
    h_d_mu2p_tt_vs_dt_mu2p->GetXaxis()->SetTitle("#Delta t_{#mu2p} (s)");
    h_d_mu2p_tt_vs_dt_mu2p->GetXaxis()->CenterTitle(kTRUE);
    h_d_mu2p_tt_vs_dt_mu2p->GetYaxis()->SetTitle("d_{#mu2p} (cm)");
    h_d_mu2p_tt_vs_dt_mu2p->GetYaxis()->CenterTitle(kTRUE);
    h_d_mu2p_tt_vs_dt_mu2p->GetYaxis()->SetTitleOffset(1.5);
    h_d_mu2p_tt_vs_dt_mu2p->Draw("COLZ");
    c_d_mu2p_tt_vs_dt_mu2p->SetTickx();
    c_d_mu2p_tt_vs_dt_mu2p->SetTicky();
    c_d_mu2p_tt_vs_dt_mu2p->Update();
}

void daq_time(const std::string& filename) {
    TChain* chain = new TChain("DAQTree");
    if (!chain) {
        std::cerr << "Cannot create TChain DAQTree\n";
        return;
    }
    chain->Add(filename.c_str());
    time_t daq_sec;
    int daq_nsec;
    time_t muveto_sec;
    int muveto_nsec;
    chain->SetBranchAddress("daq_sec", &daq_sec);
    chain->SetBranchAddress("daq_nsec", &daq_nsec);
    chain->SetBranchAddress("muveto_sec", &muveto_sec);
    chain->SetBranchAddress("muveto_nsec", &muveto_nsec);
    timestamp tot_ts, tot_ts_mu;
    for (int k = 0; k < chain->GetEntries(); ++k) {
        chain->GetEntry(k);
        timestamp ts{daq_sec, daq_nsec};
        timestamp ts_mu{muveto_sec, muveto_nsec};
        tot_ts += ts;
        tot_ts_mu += ts_mu;
    }
    std::cout << "Total DAQ time: " << tot_ts << '\n';
    std::cout << "Total MuVeto time: " << tot_ts_mu << '\n';
    double tot_seconds = timestamp_to_double(tot_ts);
    double tot_seconds_mu = timestamp_to_double(tot_ts_mu);
    std::cout << "Total DAQ time in days: " << tot_seconds / (3600.0 * 24.0) << '\n';
    std::cout << "Total MuVeto time in days: " << tot_seconds_mu / (3600.0 * 24.0) << '\n';
}

void print_all_entries(const std::string& filename, analysis_base* analysis) {
    if (!analysis->retrieve(filename)) return;
    std::set<IBD> ibds;
    std::cout << "=== Analysis: " << analysis->name << " (Total Entries: " << analysis->size() << ") ===\n";
    for (long k = 0; k < analysis->size(); ++k) {
        analysis->entry(k);
        if (k % 1000 == 0) {
            std::cout << "\rProcessing Entry " << k << " / " << analysis->size();
        }
        // analysis->print();
        if (!analysis->selection()) continue;
        ibds.insert({analysis->prompt, analysis->delayed});
    }

    for (std::set<IBD>::const_iterator it = ibds.begin(); it != ibds.end(); ++it) {
        std::cout << "Prompt: " << "E = " << it->prompt.e << ", Q = " << it->prompt.q << ", Time = " << it->prompt.ts << '\n';
        std::cout << "Delayed: " << "E = " << it->delayed.e << ", Q = " << it->delayed.q << ", Time = " << it->delayed.ts << '\n';
    }
}

struct VanessaIBD {

    int run_id;

    timestamp ts_p;
    double e_p;
    double totq_p;

    timestamp ts_d;
    double e_d;
    double totq_d;

};

bool operator<(const VanessaIBD& lhs, const VanessaIBD& rhs) {
    return lhs.ts_p < rhs.ts_p;
}

std::vector<VanessaIBD> analyze_vanessa_result() {
    TFile* file = TFile::Open("/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/ibd/summary/ReProd25C/IBD_all_reprodC.root", "READ");
    TTree* tree = file->Get<TTree>("events");
    std::set<VanessaIBD> ibds;
    double run_id;
    double t_p, t_d;
    double e_p, e_d;
    double totq_p, totq_d;
    tree->SetBranchAddress("run_number", &run_id);
    tree->SetBranchAddress("time_p_ns", &t_p);
    tree->SetBranchAddress("time_d_ns", &t_d);
    tree->SetBranchAddress("energy_p_omilrec", &e_p);
    tree->SetBranchAddress("energy_d_omilrec", &e_d);
    tree->SetBranchAddress("NPE_p", &totq_p);
    tree->SetBranchAddress("NPE_d", &totq_d);
    for (long k = 0l; k < tree->GetEntries(); ++k) {
        tree->GetEntry(k);
        time_t sec_p = static_cast<time_t>(t_p / 1.0e9);
        int nsec_p = static_cast<int>(t_p - static_cast<double>(sec_p) * 1.0e9);
        time_t sec_d = static_cast<time_t>(t_d / 1.0e9);
        int nsec_d = static_cast<int>(t_d - static_cast<double>(sec_d) * 1.0e9);
        VanessaIBD ibd;
        ibd.run_id = run_id;
        ibd.ts_p = timestamp{sec_p, nsec_p};
        ibd.e_p = e_p;
        ibd.totq_p = totq_p;
        ibd.ts_d = timestamp{sec_d, nsec_d};
        ibd.e_d = e_d;
        ibd.totq_d = totq_d;
        ibds.insert(ibd);
    }
    std::vector<VanessaIBD> ibds_vector;
    ibds_vector.reserve(ibds.size());
    for (std::set<VanessaIBD>::const_iterator it = ibds.begin(); it != ibds.end(); ++it) {
        ibds_vector.push_back(*it);
    }
    return ibds_vector;
}

void compare_with_vanessa(const std::string& filename, ibd_analysis* analysis) {
    TFile* vanessa_file = TFile::Open("/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/ibd/summary/ReProd25C/IBD_all_reprodC.root", "READ");
    TTree* vanessa_tree = vanessa_file->Get<TTree>("events");
    std::set<VanessaIBD> vanessa_ibds;
    double run_id;
    double t_p, t_d;
    double e_p, e_d;
    double totq_p, totq_d;
    vanessa_tree->SetBranchAddress("run_number", &run_id);
    vanessa_tree->SetBranchAddress("time_p_ns", &t_p);
    vanessa_tree->SetBranchAddress("time_d_ns", &t_d);
    vanessa_tree->SetBranchAddress("energy_p_omilrec", &e_p);
    vanessa_tree->SetBranchAddress("energy_d_omilrec", &e_d);
    vanessa_tree->SetBranchAddress("NPE_p", &totq_p);
    vanessa_tree->SetBranchAddress("NPE_d", &totq_d);
    for (long k = 0l; k < vanessa_tree->GetEntries(); ++k) {
        vanessa_tree->GetEntry(k);
        time_t sec_p = static_cast<time_t>(t_p / 1.0e9);
        int nsec_p = static_cast<int>(t_p - static_cast<double>(sec_p) * 1.0e9);
        time_t sec_d = static_cast<time_t>(t_d / 1.0e9);
        int nsec_d = static_cast<int>(t_d - static_cast<double>(sec_d) * 1.0e9);
        VanessaIBD ibd;
        ibd.run_id = run_id;
        ibd.ts_p = timestamp{sec_p, nsec_p};
        ibd.e_p = e_p;
        ibd.totq_p = totq_p;
        ibd.ts_d = timestamp{sec_d, nsec_d};
        ibd.e_d = e_d;
        ibd.totq_d = totq_d;
        vanessa_ibds.insert(ibd);
    }

    if (!analysis->retrieve(filename)) return;
    std::set<IBD> ibds;
    std::cout << "=== Analysis: " << analysis->name << " (Total Entries: " << analysis->size() << ") ===\n";
    for (long k = 0; k < analysis->size(); ++k) {
        analysis->entry(k);
        if (k % 1000 == 0) {
            std::cout << "\rProcessing Entry " << k << " / " << analysis->size();
        }
        std::set<VanessaIBD>::const_iterator it = std::find_if(
            vanessa_ibds.begin(),
            vanessa_ibds.end(),
            [&](const VanessaIBD& vanessa_ibd) {
                return (
                    vanessa_ibd.ts_p.sec == analysis->prompt.ts.sec &&
                    vanessa_ibd.ts_d.sec == analysis->delayed.ts.sec &&
                    vanessa_ibd.e_p == analysis->prompt.e &&
                    vanessa_ibd.e_d == analysis->delayed.e
                );
            }
        );
        bool is_only_in_vanessa = (it != vanessa_ibds.end() && !analysis->selection());
        bool is_only_in_analysis = (it == vanessa_ibds.end() && analysis->selection());
        if (!is_only_in_vanessa && !is_only_in_analysis) continue;

        std::size_t nb_neutron_veto = 0ul;
        for (std::size_t k = 0ul; k < analysis->e_n->size(); ++k) {
            if (analysis->e_n->operator[](k) < 1.5 || 20.0 < analysis->e_n->operator[](k)) continue;
            timestamp ts_n{analysis->sec_n->operator[](k), analysis->nsec_n->operator[](k)};
            vec3 pos_n{analysis->posx_n->operator[](k), analysis->posy_n->operator[](k), analysis->posz_n->operator[](k)};
            // if (pos_n.Mag() > 17700.0) continue;
            if (mag(analysis->prompt.pos - pos_n) > 4000.0 || mag(analysis->delayed.pos - pos_n) > 4000.0) continue;
            if (analysis->prompt.ts - ts_n < timestamp{0, 20000} || timestamp{0, 1200000000} < analysis->prompt.ts - ts_n) continue;
            if (analysis->delayed.ts - ts_n < timestamp{0, 20000} || timestamp{0, 1200000000} < analysis->delayed.ts - ts_n) continue;
            ++nb_neutron_veto;
        }

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < analysis->e_mult->size(); ++k) {
            if (analysis->e_mult->operator[](k) < 2.0 || 12.0 < analysis->e_mult->operator[](k)) continue;
            timestamp ts_mult{analysis->sec_mult->operator[](k), analysis->nsec_mult->operator[](k)};
            vec3 pos_mult{analysis->posx_mult->operator[](k), analysis->posy_mult->operator[](k), analysis->posz_mult->operator[](k)};
            // if (pos_mult.Mag() > 17700.0) continue;
            // if ((analysis->prompt.pos - pos_mult).Mag() > 4000.0 || (analysis->delayed.pos - pos_mult).Mag() > 4000.0) continue;
            if (ts_mult < analysis->prompt.ts - timestamp{0, 1000000} || analysis->delayed.ts + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }

        std::string method = "";
        double d_mu2p = std::numeric_limits<double>::infinity();
        double t_mu2p = std::numeric_limits<double>::infinity();
        double d_mu2d = std::numeric_limits<double>::infinity();
        double t_mu2d = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0ul; k < analysis->method_mu->size(); ++k) {
            timestamp ts_mu{analysis->sec_mu->operator[](k), analysis->nsec_mu->operator[](k)};
            vec3 pos_mu{analysis->posx_mu->operator[](k), analysis->posy_mu->operator[](k), analysis->posz_mu->operator[](k)};
            vec3 dir_mu{analysis->dirx_mu->operator[](k), analysis->diry_mu->operator[](k), analysis->dirz_mu->operator[](k)};
            double tmp_d_mu2p = mag(cross(dir_mu, analysis->prompt.pos - pos_mu));
            double tmp_t_mu2p = timestamp_to_double(analysis->prompt.ts - ts_mu);
            double tmp_d_mu2d = mag(cross(dir_mu, analysis->delayed.pos - pos_mu));
            double tmp_t_mu2d = timestamp_to_double(analysis->delayed.ts - ts_mu);
            if (tmp_d_mu2p < d_mu2p && 0.0 < tmp_t_mu2p && tmp_t_mu2p < 1.2) {
                method = analysis->method_mu->operator[](k);
                d_mu2p = tmp_d_mu2p;
                t_mu2p = tmp_t_mu2p;
                d_mu2d = tmp_d_mu2d;
                t_mu2d = tmp_t_mu2d;
            }
        }
        
        std::cout << (is_only_in_vanessa ? "[Only in Vanessa] ======================================" : "");
        std::cout << (is_only_in_analysis ? "[Only in Analysis] ======================================" : "");
        std::cout << '\n';
        analysis->print();
        std::cout << "Number of Neutron Veto associated: " << nb_neutron_veto << '\n';
        std::cout << "Number of Multiplicity Veto associated: " << nb_multu_veto << '\n';

        // print neutrons
        for (std::size_t k = 0ul; k < analysis->e_n->size(); ++k) {
            timestamp ts_n{analysis->sec_n->operator[](k), analysis->nsec_n->operator[](k)};
            vec3 pos_n{analysis->posx_n->operator[](k), analysis->posy_n->operator[](k), analysis->posz_n->operator[](k)};
            double e_n = analysis->e_n->operator[](k);
            double totq_n = analysis->totq_n->operator[](k);
            std::cout << "  Neutron: " << pos_n << ", E = " << e_n << ", Q = " << totq_n << ", Time = " << ts_n << '\n';
        }

        // print mults
        for (std::size_t k = 0ul; k < analysis->e_mult->size(); ++k) {
            timestamp ts_mult{analysis->sec_mult->operator[](k), analysis->nsec_mult->operator[](k)};
            vec3 pos_mult{analysis->posx_mult->operator[](k), analysis->posy_mult->operator[](k), analysis->posz_mult->operator[](k)};
            double e_mult = analysis->e_mult->operator[](k);
            double totq_mult = analysis->totq_mult->operator[](k);
            std::cout << "  Mult: " << pos_mult << ", E = " << e_mult << ", Q = " << totq_mult << ", Time = " << ts_mult << '\n';
        }

        // print muons
        for (std::size_t k = 0ul; k < analysis->method_mu->size(); ++k) {
            timestamp ts_mu{analysis->sec_mu->operator[](k), analysis->nsec_mu->operator[](k)};
            vec3 pos_mu{analysis->posx_mu->operator[](k), analysis->posy_mu->operator[](k), analysis->posz_mu->operator[](k)};
            vec3 dir_mu{analysis->dirx_mu->operator[](k), analysis->diry_mu->operator[](k), analysis->dirz_mu->operator[](k)};
            double totq_mu = analysis->totq_mu->operator[](k);
            double quality_mu = analysis->quality_mu->operator[](k);
            std::cout << "  Muon (" << analysis->method_mu->operator[](k) << "): " << pos_mu << ", Dir = " << dir_mu << ", Q = " << totq_mu << ", Quality = " << quality_mu << ", Time = " << ts_mu << '\n';
        }
    }
}

std::vector<double> linspace_cpp(double start, double stop, int num) {
    if (num <= 1) {
        if (num == 1) return {start};
        std::cerr << "Warning: linspace_cpp requires num >= 1. Returning empty vector\n";
        return {};
    }
    double step = (stop - start) / (num - 1);
    std::vector<double> result;
    for (int i = 0; i < num; ++i) {
        double value = start + static_cast<double>(i) * step;
        if (i == num - 1) {
            result.push_back(stop);
        } 
        else {
            result.push_back(value);
        }
    }
    
    return result;
}

std::vector<double> generate_segment_boundaries(double start, double stop, int num_bins) {
    if (num_bins <= 0) return {};
    int num_points = num_bins + 1;
    double expected_width = (stop - start) / num_bins;
    
    std::vector<double> segment;
    segment.reserve(num_points);
    segment.push_back(start); 

    for (int i = 1; i < num_points; ++i) {
        double boundary = start + i * expected_width;
        
        if (i == num_points - 1) {
             segment.push_back(stop);
        } else {
             segment.push_back(boundary);
        }
    }
    return segment;
}

std::vector<double> create_custom_e_p_bins() {

    // double s1_start = 0.8;
    // double s2_start = 0.94;
    // double s3_start = 7.44;
    // double s4_start = 7.8;
    // double s5_start = 8.2;
    // double stop = 12.0;

    // int s1_bins = 1;
    // int s2_bins = 325;
    // int s3_bins = 9;
    // int s4_bins = 4;
    // int s5_bins = 1;
    // int tot_bins = s1_bins + s2_bins + s3_bins + s4_bins + s5_bins;

    double edges[] = {0.7, 1.0, 6.6, 7.4, 7.7, 8.1, 8.6, 9.4, 12.0};
    int    bins[]  = {  1,  56,   4,   1,   1,   1,   1,   1};
    
    std::vector<double> e_p_bins;

    for (std::size_t k = 0ul; k < 8ul; ++k) {
        double start = edges[k];
        double stop = edges[k + 1];
        int nbins = bins[k];
        std::vector<double> segment = generate_segment_boundaries(start, stop, nbins);
        if (k == 0ul) {
            e_p_bins.insert(e_p_bins.end(), segment.begin(), segment.end());
        } else {
            e_p_bins.insert(e_p_bins.end(), segment.begin() + 1, segment.end());
        }
    }
    
    return e_p_bins;
}

#undef PRINT_ALL_ENTRIES
#undef COMPARE_WITH_VANESSA
#define GET_ALL_IBD

void analysisgroupc_printer(const std::string& filename, const std::string& suffix) {
    daq_time(filename);

    cosmo_rate_analysis cosmo_rate_with_neutron_analysis(suffix);
    analyze_cosmo_rate_with_neutron(filename, &cosmo_rate_with_neutron_analysis);

    cosmo_shape_analysis cosmo_before_analysis(suffix, timestamp{0, -1200000000}, timestamp{0, -5000000}, 3000.0);
    cosmo_shape_analysis cosmo_after_analysis(suffix, timestamp{0, 5000000}, timestamp{0, 1200000000}, 3000.0);
    std::vector<Cosmo> cosmos_before = get_all_cosmo(filename, &cosmo_before_analysis);
    std::vector<Cosmo> cosmos_after = get_all_cosmo(filename, &cosmo_after_analysis);

    ibd_analysis ibd_analysis(suffix);
#ifdef PRINT_ALL_ENTRIES
    print_all_entries(filename, &ibd_analysis);
#endif 
#ifdef COMPARE_WITH_VANESSA
    compare_with_vanessa(filename, &ibd_analysis);
#endif
#ifdef GET_ALL_IBD
    std::vector<IBD> ibds = get_all_ibd(filename, &ibd_analysis);
    std::vector<VanessaIBD> vanessa_ibds = analyze_vanessa_result();


    // double e_p_min = 0.7;
    // double e_p_max = 12.0;
    // double e_p_width = 0.20;
    // int e_p_nbin = std::round((e_p_max - e_p_min) / e_p_width) + 1;
    // std::vector<double> e_p_bins = linspace_cpp(e_p_min, e_p_max, e_p_nbin);
    std::vector<double> e_p_bins = create_custom_e_p_bins();
    TH1D* h_e_p = new TH1D("h_e_p", "Prompt energy", e_p_bins.size() - 1, e_p_bins.data());
    TH1D* h_e_p_vanessa = new TH1D("h_e_p_vanessa", "Prompt energy (Vanessa)", e_p_bins.size() - 1, e_p_bins.data());
    TH1D* h_e_p_cosmo_before = new TH1D("h_e_p_cosmo_before", "Prompt energy (Cosmo Before)", e_p_bins.size() - 1, e_p_bins.data());
    TH1D* h_e_p_cosmo_after = new TH1D("h_e_p_cosmo_after", "Prompt energy (Cosmo After)", e_p_bins.size() - 1, e_p_bins.data());
    TH1D* h_e_p_cosmo_diff = new TH1D("h_e_p_cosmo_diff", "Prompt energy (Cosmo After - Cosmo Before)", e_p_bins.size() - 1, e_p_bins.data());

    double e_d_min = 2.0;
    double e_d_max = 2.5;
    double e_d_width = 0.02;
    int e_d_nbin = std::round((e_d_max - e_d_min) / e_d_width) + 1;
    std::vector<double> e_d_bins = linspace_cpp(e_d_min, e_d_max, e_d_nbin);
    TH1D* h_e_d = new TH1D("h_e_d", "Delayed energy", e_d_bins.size() - 1, e_d_bins.data());
    TH1D* h_e_d_vanessa = new TH1D("h_e_d_vanessa", "Delayed energy (Vanessa)", e_d_bins.size() - 1, e_d_bins.data());
    TH1D* h_e_d_cosmo_before = new TH1D("h_e_d_cosmo_before", "Delayed energy (Cosmo Before)", e_d_bins.size() - 1, e_d_bins.data());
    TH1D* h_e_d_cosmo_after = new TH1D("h_e_d_cosmo_after", "Delayed energy (Cosmo After)", e_d_bins.size() - 1, e_d_bins.data());

    double e_dt_min = 0.0;
    double e_dt_max = 1.0;
    double e_dt_width = 0.025;
    int e_dt_nbin = std::round((e_dt_max - e_dt_min) / e_dt_width) + 1;
    std::vector<double> e_dt_bins = linspace_cpp(e_dt_min, e_dt_max, e_dt_nbin);
    TH1D* h_dt = new TH1D("h_dt", "Prompt-Delayed time difference", e_dt_bins.size() - 1, e_dt_bins.data());
    TH1D* h_dt_vanessa = new TH1D("h_dt_vanessa", "Prompt-Delayed time difference (Vanessa)", e_dt_bins.size() - 1, e_dt_bins.data());
    TH1D* h_dt_cosmo_before = new TH1D("h_dt_cosmo_before", "Prompt-Delayed time difference (Cosmo Before)", e_dt_bins.size() - 1, e_dt_bins.data());
    TH1D* h_dt_cosmo_after = new TH1D("h_dt_cosmo_after", "Prompt-Delayed time difference (Cosmo After)", e_dt_bins.size() - 1, e_dt_bins.data());

    double e_dr_min = 0.0;
    double e_dr_max = 1.5;
    double e_dr_width = 0.05;
    int e_dr_nbin = std::round((e_dr_max - e_dr_min) / e_dr_width) + 1;
    std::vector<double> e_dr_bins = linspace_cpp(e_dr_min, e_dr_max, e_dr_nbin);
    TH1D* h_dr = new TH1D("h_dr", "Prompt-Delayed distance", e_dr_bins.size() - 1, e_dr_bins.data());
    TH1D* h_dr_cosmo_before = new TH1D("h_dr_cosmo_before", "Prompt-Delayed distance (Cosmo Before)", e_dr_bins.size() - 1, e_dr_bins.data());
    TH1D* h_dr_cosmo_after = new TH1D("h_dr_cosmo_after", "Prompt-Delayed distance (Cosmo After)", e_dr_bins.size() - 1, e_dr_bins.data());

    double rho_min = 0.0;
    double rho_max = 17.7 * 17.7;
    int rho_nbin = 51;
    double z_min = -20.0;
    double z_max = 20.0;
    int z_nbin = 51;
    std::vector<double> rho_bins = linspace_cpp(rho_min, rho_max, rho_nbin);
    std::vector<double> z_bins = linspace_cpp(z_min, z_max, z_nbin);
    TH2D* h_rho_z_p = new TH2D("h_rho_z_p", "Prompt vertex distribution", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
    TH2D* h_rho_z_d = new TH2D("h_rho_z_d", "Delayed vertex distribution", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
    TH2D* h_rho_z_p_cosmo_before = new TH2D("h_rho_z_p_cosmo_before", "Prompt vertex distribution (Cosmo Before)", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
    TH2D* h_rho_z_d_cosmo_before = new TH2D("h_rho_z_d_cosmo_before", "Delayed vertex distribution (Cosmo Before)", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
    TH2D* h_rho_z_p_cosmo_after = new TH2D("h_rho_z_p_cosmo_after", "Prompt vertex distribution (Cosmo After)", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
    TH2D* h_rho_z_d_cosmo_after = new TH2D("h_rho_z_d_cosmo_after", "Delayed vertex distribution (Cosmo After)", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());

    for (const IBD& ibd : ibds) {
        h_e_p->Fill(ibd.prompt.e);
        h_e_d->Fill(ibd.delayed.e);
        h_dt->Fill(timestamp_to_double(ibd.delayed.ts - ibd.prompt.ts) * 1000.0);
        h_dr->Fill(mag(ibd.delayed.pos - ibd.prompt.pos) / 1000.0);
        h_rho_z_p->Fill((ibd.prompt.pos.x * ibd.prompt.pos.x + ibd.prompt.pos.y * ibd.prompt.pos.y) / 1.0e6, ibd.prompt.pos.z / 1000.0);
        h_rho_z_d->Fill((ibd.delayed.pos.x * ibd.delayed.pos.x + ibd.delayed.pos.y * ibd.delayed.pos.y) / 1.0e6, ibd.delayed.pos.z / 1000.0);
    }

    for (const Cosmo& cosmo : cosmos_before) {
        h_e_p_cosmo_before->Fill(cosmo.prompt.e);
        h_e_d_cosmo_before->Fill(cosmo.delayed.e);
        h_dt_cosmo_before->Fill(timestamp_to_double(cosmo.delayed.ts - cosmo.prompt.ts) * 1000.0);
        h_dr_cosmo_before->Fill(mag(cosmo.delayed.pos - cosmo.prompt.pos) / 1000.0);
        h_rho_z_p_cosmo_before->Fill((cosmo.prompt.pos.x * cosmo.prompt.pos.x + cosmo.prompt.pos.y * cosmo.prompt.pos.y) / 1.0e6, cosmo.prompt.pos.z / 1000.0);
        h_rho_z_d_cosmo_before->Fill((cosmo.delayed.pos.x * cosmo.delayed.pos.x + cosmo.delayed.pos.y * cosmo.delayed.pos.y) / 1.0e6, cosmo.delayed.pos.z / 1000.0);
    }

    for (const Cosmo& cosmo : cosmos_after) {
        h_e_p_cosmo_after->Fill(cosmo.prompt.e);
        h_e_d_cosmo_after->Fill(cosmo.delayed.e);
        h_dt_cosmo_after->Fill(timestamp_to_double(cosmo.delayed.ts - cosmo.prompt.ts) * 1000.0);
        h_dr_cosmo_after->Fill(mag(cosmo.delayed.pos - cosmo.prompt.pos) / 1000.0);
        h_rho_z_p_cosmo_after->Fill((cosmo.prompt.pos.x * cosmo.prompt.pos.x + cosmo.prompt.pos.y * cosmo.prompt.pos.y) / 1.0e6, cosmo.prompt.pos.z / 1000.0);
        h_rho_z_d_cosmo_after->Fill((cosmo.delayed.pos.x * cosmo.delayed.pos.x + cosmo.delayed.pos.y * cosmo.delayed.pos.y) / 1.0e6, cosmo.delayed.pos.z / 1000.0);
    }

    h_e_p_cosmo_diff->Add(h_e_p_cosmo_after, h_e_p_cosmo_before, 1.0, -1.0);

    for (const VanessaIBD& ibd : vanessa_ibds) {
        h_e_p_vanessa->Fill(ibd.e_p);
        h_e_d_vanessa->Fill(ibd.e_d);
        h_dt_vanessa->Fill(timestamp_to_double(ibd.ts_d - ibd.ts_p) * 1000.0);
        // h_dr_vanessa->Fill((ibd.delayed.pos - ibd.prompt.pos).Mag() / 1000.0);
        // h_rho_z_p_vanessa->Fill((ibd.prompt.pos.X() * ibd.prompt.pos.X() + ibd.prompt.pos.Y() * ibd.prompt.pos.Y()) / 1.0e6, ibd.prompt.pos.Z() / 1000.0);
        // h_rho_z_d_vanessa->Fill((ibd.delayed.pos.X() * ibd.delayed.pos.X() + ibd.delayed.pos.Y() * ibd.delayed.pos.Y()) / 1.0e6, ibd.delayed.pos.Z() / 1000.0);
    }

    // ============================================================================================
    // Prompt energy
    // ============================================================================================

    TCanvas* c_e_p = new TCanvas("c_e_p", "Prompt energy", 1000, 1000);
    c_e_p->cd();

    h_e_p->SetLineWidth(3);
    h_e_p->SetLineStyle(kSolid);
    h_e_p->SetLineColorAlpha(kBlue, 1.0);

    h_e_p_vanessa->SetLineWidth(3);
    h_e_p_vanessa->SetLineStyle(kSolid);
    h_e_p_vanessa->SetLineColorAlpha(kRed, 1.0);

    h_e_p->Draw();
    h_e_p_vanessa->Draw("SAME");

    c_e_p->Update();

    // ============================================================================================
    // Delayed energy
    // ============================================================================================

    TCanvas* c_e_d = new TCanvas("c_e_d", "Delayed energy", 1000, 1000);
    c_e_d->cd();

    h_e_d->SetLineWidth(3);
    h_e_d->SetLineStyle(kSolid);
    h_e_d->SetLineColorAlpha(kBlue, 1.0);

    h_e_d_vanessa->SetLineWidth(3);
    h_e_d_vanessa->SetLineStyle(kSolid);
    h_e_d_vanessa->SetLineColorAlpha(kRed, 1.0);

    h_e_d->Draw();
    h_e_d_vanessa->Draw("SAME");

    c_e_d->Update();

    // ============================================================================================
    // Prompt-Delayed time difference
    // ============================================================================================

    TCanvas* c_dt = new TCanvas("c_dt", "Prompt-Delayed time difference", 1000, 1000);
    c_dt->cd();

    h_dt->SetLineWidth(3);
    h_dt->SetLineStyle(kSolid);
    h_dt->SetLineColorAlpha(kBlue, 1.0);

    h_dt_vanessa->SetLineWidth(3);
    h_dt_vanessa->SetLineStyle(kSolid);
    h_dt_vanessa->SetLineColorAlpha(kRed, 1.0);

    h_dt->Draw();
    h_dt_vanessa->Draw("SAME");
    
    c_dt->Update();

    // ============================================================================================
    // Prompt-Delayed distance
    // ============================================================================================

    TCanvas* c_dr = new TCanvas("c_dr", "Prompt-Delayed distance", 1000, 1000);
    c_dr->cd();
    h_dr->Draw();
    c_dr->Update();

    // ============================================================================================
    // Prompt vertex position
    // ============================================================================================

    TCanvas* c_rho_z_p = new TCanvas("c_rho_z_p", "Prompt vertex distribution", 1000, 1000);
    c_rho_z_p->cd();
    h_rho_z_p->Draw();
    c_rho_z_p->Update();

    // ============================================================================================
    // Delayed vertex position
    // ============================================================================================

    TCanvas* c_rho_z_d = new TCanvas("c_rho_z_d", "Delayed vertex distribution", 1000, 1000);
    c_rho_z_d->cd();
    h_rho_z_d->Draw();
    c_rho_z_d->Update();
#endif

    // ============================================================================================
    // Cosmo before - Prompt energy
    // ============================================================================================

    TCanvas* c_e_p_cosmo_before = new TCanvas("c_e_p_cosmo_before", "Prompt energy (Cosmo before)", 1000, 1000);
    c_e_p_cosmo_before->cd();

    h_e_p_cosmo_before->SetLineWidth(3);
    h_e_p_cosmo_before->SetLineStyle(kSolid);
    h_e_p_cosmo_before->SetLineColorAlpha(kBlue, 1.0);

    h_e_p_cosmo_before->Draw();

    c_e_p_cosmo_before->Update();

    // ============================================================================================
    // Cosmo before - Delayed energy
    // ============================================================================================

    TCanvas* c_e_d_cosmo_before = new TCanvas("c_e_d_cosmo_before", "Delayed energy (Cosmo before)", 1000, 1000);
    c_e_d_cosmo_before->cd();

    h_e_d_cosmo_before->SetLineWidth(3);
    h_e_d_cosmo_before->SetLineStyle(kSolid);
    h_e_d_cosmo_before->SetLineColorAlpha(kBlue, 1.0);

    h_e_d_cosmo_before->Draw();

    c_e_d_cosmo_before->Update();

    // ============================================================================================
    // Cosmo before - Prompt-Delayed time difference
    // ============================================================================================

    TCanvas* c_dt_cosmo_before = new TCanvas("c_dt_cosmo_before", "Prompt-Delayed time difference (Cosmo before)", 1000, 1000);
    c_dt_cosmo_before->cd();
    
    h_dt_cosmo_before->SetLineWidth(3);
    h_dt_cosmo_before->SetLineStyle(kSolid);
    h_dt_cosmo_before->SetLineColorAlpha(kBlue, 1.0);

    h_dt_cosmo_before->Draw();

    c_dt_cosmo_before->Update();

    // ============================================================================================
    // Cosmo before - Prompt-Delayed distance
    // ============================================================================================

    TCanvas* c_dr_cosmo_before = new TCanvas("c_dr_cosmo_before", "Prompt-Delayed distance (Cosmo before)", 1000, 1000);
    c_dr_cosmo_before->cd();

    h_dr_cosmo_before->Draw();

    c_dr_cosmo_before->Update();

    // ============================================================================================
    // Cosmo before - Prompt vertex position
    // ============================================================================================

    TCanvas* c_rho_z_p_cosmo_before = new TCanvas("c_rho_z_p_cosmo_before", "Prompt vertex distribution (Cosmo before)", 1000, 1000);
    c_rho_z_p_cosmo_before->cd();

    h_rho_z_p_cosmo_before->Draw();

    c_rho_z_p_cosmo_before->Update();

    // ============================================================================================
    // Cosmo before - Delayed vertex position
    // ============================================================================================

    TCanvas* c_rho_z_d_cosmo_before = new TCanvas("c_rho_z_d_cosmo_before", "Delayed vertex distribution (Cosmo before)", 1000, 1000);
    c_rho_z_d_cosmo_before->cd();

    h_rho_z_d_cosmo_before->Draw();

    c_rho_z_d_cosmo_before->Update();

    // ============================================================================================
    // Cosmo after - Prompt energy
    // ============================================================================================

    TCanvas* c_e_p_cosmo_after = new TCanvas("c_e_p_cosmo_after", "Prompt energy (Cosmo after)", 1000, 1000);
    c_e_p_cosmo_after->cd();

    h_e_p_cosmo_after->SetLineWidth(3);
    h_e_p_cosmo_after->SetLineStyle(kSolid);
    h_e_p_cosmo_after->SetLineColorAlpha(kBlue, 1.0);

    h_e_p_cosmo_after->Draw();

    c_e_p_cosmo_after->Update();

    // ============================================================================================
    // Cosmo after - Delayed energy
    // ============================================================================================

    TCanvas* c_e_d_cosmo_after = new TCanvas("c_e_d_cosmo_after", "Delayed energy (Cosmo after)", 1000, 1000);
    c_e_d_cosmo_after->cd();
    
    h_e_d_cosmo_after->SetLineWidth(3);
    h_e_d_cosmo_after->SetLineStyle(kSolid);
    h_e_d_cosmo_after->SetLineColorAlpha(kBlue, 1.0);

    h_e_d_cosmo_after->Draw();

    c_e_d_cosmo_after->Update();

    // ============================================================================================
    // Cosmo after - Prompt-Delayed time difference
    // ============================================================================================

    TCanvas* c_dt_cosmo_after = new TCanvas("c_dt_cosmo_after", "Prompt-Delayed time difference (Cosmo after)", 1000, 1000);
    c_dt_cosmo_after->cd();

    h_dt_cosmo_after->SetLineWidth(3);
    h_dt_cosmo_after->SetLineStyle(kSolid);
    h_dt_cosmo_after->SetLineColorAlpha(kBlue, 1.0);

    h_dt_cosmo_after->Draw();

    c_dt_cosmo_after->Update();

    // ============================================================================================
    // Cosmo after - Prompt-Delayed distance
    // ============================================================================================

    TCanvas* c_dr_cosmo_after = new TCanvas("c_dr_cosmo_after", "Prompt-Delayed distance (Cosmo after)", 1000, 1000);
    c_dr_cosmo_after->cd();

    h_dr_cosmo_after->Draw();

    c_dr_cosmo_after->Update();

    // ============================================================================================
    // Cosmo after - Prompt vertex position
    // ============================================================================================

    TCanvas* c_rho_z_p_cosmo_after = new TCanvas("c_rho_z_p_cosmo_after", "Prompt vertex distribution (Cosmo after)", 1000, 1000);
    c_rho_z_p_cosmo_after->cd();

    h_rho_z_p_cosmo_after->Draw();

    c_rho_z_p_cosmo_after->Update();

    // ============================================================================================
    // Cosmo after - Delayed vertex position
    // ============================================================================================

    TCanvas* c_rho_z_d_cosmo_after = new TCanvas("c_rho_z_d_cosmo_after", "Delayed vertex distribution (Cosmo after)", 1000, 1000);
    c_rho_z_d_cosmo_after->cd();

    h_rho_z_d_cosmo_after->Draw();

    c_rho_z_d_cosmo_after->Update();

    // ============================================================================================
    // Cosmo diff - Prompt energy
    // ============================================================================================

    TCanvas* c_e_p_cosmo_diff = new TCanvas("c_e_p_cosmo_diff", "Prompt energy (Cosmo diff)", 1000, 1000);
    c_e_p_cosmo_diff->cd();

    h_e_p_cosmo_diff->SetLineWidth(3);
    h_e_p_cosmo_diff->SetLineStyle(kSolid);
    h_e_p_cosmo_diff->SetLineColorAlpha(kBlue, 1.0);

    h_e_p_cosmo_diff->Draw();

    c_e_p_cosmo_diff->Update();

}