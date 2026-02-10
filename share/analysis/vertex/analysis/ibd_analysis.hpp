#ifndef ANALYSIS_IBD_ANALYSIS_HPP_
#define ANALYSIS_IBD_ANALYSIS_HPP_

#include <set>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>

#include "analysis/basic_analysis.hpp"
#include "utils/plot.hpp"

class ibd_analysis : public basic_analysis {

public:

    ibd_analysis(const std::string& name, const std::string& filepath, const std::string& suffix) :
        basic_analysis{name, filepath, suffix}
    {}

    ~ibd_analysis() override = default;

    bool selection() override {
        double e_p = m_nav->prompt.e * m_gtc.interpolate(m_nav->prompt.ts);
        double e_d = m_nav->delayed.e * m_gtc.interpolate(m_nav->delayed.ts);

        if (e_p < 0.7 || 12.0 < e_p) return false;
        if (e_d < 2.0 || 2.5 < e_d) return false;
        if (mag(m_nav->prompt.pos) > 16500.0) return false;
        if (std::abs(m_nav->prompt.pos.z) > 15500.0 && std::sqrt(m_nav->prompt.pos.x * m_nav->prompt.pos.x + m_nav->prompt.pos.y * m_nav->prompt.pos.y) < 2000.0) return false;
        timestamp ts_diff = m_nav->delayed.ts - m_nav->prompt.ts;
        if (ts_diff < timestamp{0, 5000} || timestamp{0, 1000000} < ts_diff) return false;
        vec3 pos_diff = m_nav->delayed.pos - m_nav->prompt.pos;
        if (mag(pos_diff) > 1500.0) return false;

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < m_nav->e_mult.size(); ++k) {
            timestamp ts_mult{m_nav->sec_mult[k], m_nav->nsec_mult[k]};
            vec3 pos_mult{m_nav->posx_mult[k], m_nav->posy_mult[k], m_nav->posz_mult[k]};
            double e_mult = m_nav->e_mult[k] * m_gtc.interpolate(ts_mult);
            if (e_mult < 2.0 || 12.0 < e_mult) continue;
            if (ts_mult < m_nav->prompt.ts - timestamp{0, 1000000} || m_nav->delayed.ts + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }
        if (nb_multu_veto) return false;

        std::size_t nb_neutron_veto = 0ul;
        for (std::size_t k = 0ul; k < m_nav->e_n.size(); ++k) {
            timestamp ts_n{m_nav->sec_n[k], m_nav->nsec_n[k]};
            vec3 pos_n{m_nav->posx_n[k], m_nav->posy_n[k], m_nav->posz_n[k]};
            double e_n = m_nav->e_n[k] * m_gtc.interpolate(ts_n);
            if (e_n < 1.5 || 20.0 < e_n) continue;
            if (m_nav->stdt_n[k] > 275.0) continue;
            if (mag(m_nav->prompt.pos - pos_n) > 4000.0) continue;
            if (m_nav->prompt.ts < ts_n + timestamp{0, 20000} || ts_n + timestamp{0, 1200000000} < m_nav->prompt.ts) continue;
            if (mag(m_nav->delayed.pos - pos_n) > 4000.0) continue;
            if (m_nav->delayed.ts < ts_n + timestamp{0, 20000} || ts_n + timestamp{0, 1200000000} < m_nav->delayed.ts) continue;
            ++nb_neutron_veto;
        }
        if (nb_neutron_veto) return false;

        if ( std::pow((m_nav->meta_prompt.stdhit - 0.55) / 0.45, 2.0) + std::pow((m_nav->meta_prompt.stdt - 170.0) / 80.0, 2.0) > 1.0 ) return false;

        return true;
    }

    bool process() override {
        m_ibds.insert({m_nav->run_id, m_nav->prompt, m_nav->delayed});
        return true;
    }

    bool save() override {
        TFile* f = TFile::Open(Form("%s.root", m_name.c_str()), "RECREATE");
        if (!f) {
            std::cerr << "Cannot open file " << m_name << ".root for writing\n";
            return false;
        }
        TTree* t = new TTree("events", "Events");
        if (!t) {
            std::cerr << "Cannot create tree events\n";
            return false;
        }
        int run_id;
        vec3 pos_p, pos_d;
        timestamp ts_p, ts_d;
        double e_p, e_d;
        t->Branch("run_id", &run_id);
        t->Branch("posx_p", &pos_p.x);
        t->Branch("posy_p", &pos_p.y);
        t->Branch("posz_p", &pos_p.z);
        t->Branch("sec_p", &ts_p.sec);
        t->Branch("nsec_p", &ts_p.nsec);
        t->Branch("e_p", &e_p);
        t->Branch("posx_d", &pos_d.x);
        t->Branch("posy_d", &pos_d.y);
        t->Branch("posz_d", &pos_d.z);
        t->Branch("sec_d", &ts_d.sec);
        t->Branch("nsec_d", &ts_d.nsec);
        t->Branch("e_d", &e_d);

        for (std::set<ibd>::const_iterator it = m_ibds.begin(); it != m_ibds.end(); ++it) {
            run_id = it->run_id;
            pos_p = it->prompt.pos;
            pos_d = it->delayed.pos;
            ts_p = it->prompt.ts;
            ts_d = it->delayed.ts;
            e_p = it->prompt.e * m_gtc.interpolate(it->prompt.ts);
            e_d = it->delayed.e * m_gtc.interpolate(it->delayed.ts);
            t->Fill();
        }

        f->cd();
        t->Write();
        f->Close();
        return true;
    }

    void result() override {
        std::vector<ibd> ibds;
        ibds.reserve(m_ibds.size());
        for (std::set<ibd>::const_iterator it = m_ibds.begin(); it != m_ibds.end(); ++it) {
            ibds.push_back(*it);
        }
    
        TH1D* h_e_p = make_prompt_energy_plot("h_e_p", "Prompt energy", ibds);
        TH1D* h_e_d = make_delayed_energy_plot("h_e_d", "Delayed energy", ibds);
        TH1D* h_dt = make_prompt_delayed_time_plot("h_dt", "Prompt-Delayed time difference", ibds);
        TH1D* h_dr = make_prompt_delayed_distance_plot("h_dr", "Prompt-Delayed distance", ibds);
        TH2D* h_rho_z_p = make_prompt_spatial_plot("h_rho_z_p", "Prompt vertex distribution", ibds);
        TH2D* h_rho_z_d = make_delayed_spatial_plot("h_rho_z_d", "Delayed vertex distribution", ibds);

        // ============================================================================================
        // Prompt energy
        // ============================================================================================

        h_e_p->SetStats(false);
        pimp_my_line(h_e_p, LineConfig{.style = kSolid, .width = 3, .color = kBlue});
        pimp_my_axis(h_e_p->GetXaxis(), AxisConfig{.ndivisions = 406, .maxdigits = 3});
        pimp_my_name(h_e_p->GetXaxis(), NameConfig{.title = "E_{d} (MeV)"});
        h_e_p->GetXaxis()->CenterTitle(true);
        pimp_my_axis(h_e_p->GetYaxis(), AxisConfig{.maxdigits = 3, .title = {.offset = 1.25}});
        pimp_my_name(h_e_p->GetYaxis(), NameConfig{.title = "Entries"});
        h_e_p->GetYaxis()->CenterTitle(true);
        TCanvas* c_e_p = plot_basic(h_e_p, "HIST");

        // ============================================================================================
        // Delayed energy
        // ============================================================================================
        
        h_e_d->SetStats(false);
        pimp_my_line(h_e_d, LineConfig{.style = kSolid, .width = 3, .color = kBlue});
        pimp_my_axis(h_e_d->GetXaxis(), AxisConfig{.ndivisions = 405, .maxdigits = 3});
        pimp_my_name(h_e_d->GetXaxis(), NameConfig{.title = "E_{d} (MeV)"});
        h_e_d->GetXaxis()->CenterTitle(true);
        pimp_my_axis(h_e_d->GetYaxis(), AxisConfig{.maxdigits = 3, .title = {.offset = 1.25}});
        pimp_my_name(h_e_d->GetYaxis(), NameConfig{.title = "Entries"});
        h_e_d->GetYaxis()->CenterTitle(true);
        TF1* f_e_d = create_gaussian_function(h_e_d, 2.0, 2.5);
        TFitResultPtr fres_e_d = h_e_d->Fit(f_e_d, "R");
        TPaveText* pt_e_d = new TPaveText(0.5, 0.65, 0.85, 0.85, "NDC");
        const double A_gaus    = f_e_d->GetParameter(0);
        const double Aerr_gaus = f_e_d->GetParError(0);
        const double mu        = f_e_d->GetParameter(1);
        const double muerr     = f_e_d->GetParError(1);
        const double sigma     = f_e_d->GetParameter(2);
        const double sigmaerr  = f_e_d->GetParError(2);
        pt_e_d->SetFillStyle(0);
        pt_e_d->SetBorderSize(0);
        pt_e_d->SetTextAlign(12);
        pt_e_d->SetTextFont(42);
        pt_e_d->SetTextSize(0.035);
        pt_e_d->AddText("#bf{Fit:} A exp(-(x - #mu)^2 / 2 #sigma^2)");
        pt_e_d->AddText(Form("A = %.3g #pm %.3g", A_gaus, Aerr_gaus));
        pt_e_d->AddText(Form("#mu = %.3g #pm %.3g MeV", mu, muerr));
        pt_e_d->AddText(Form("#sigma = %.3g #pm %.3g MeV", sigma, sigmaerr));
        TCanvas* c_e_d = plot_basic(h_e_d, "HIST");
        f_e_d->Draw("SAME");
        pt_e_d->Draw("SAME");
        c_e_p->Update();

        // ============================================================================================
        // Prompt-Delayed time difference
        // ============================================================================================

        h_dt->SetStats(false);
        pimp_my_line(h_dt, LineConfig{.style = kSolid, .width = 3, .color = kBlue});
        pimp_my_axis(h_dt->GetXaxis(), AxisConfig{.ndivisions = 405, .maxdigits = 3});
        pimp_my_name(h_dt->GetXaxis(), NameConfig{.title = "#Delta t (ms)"});
        h_dt->GetXaxis()->CenterTitle(true);
        pimp_my_axis(h_dt->GetYaxis(), AxisConfig{.maxdigits = 3, .title = {.offset = 1.25}});
        pimp_my_name(h_dt->GetYaxis(), NameConfig{.title = "Entries"});
        h_dt->GetYaxis()->CenterTitle(true);
        TF1* f_dt = create_exponential_decay_function(h_dt, 0.05, 1.0);
        TFitResultPtr fres_dt = h_dt->Fit(f_dt, "R");
        TPaveText* pt_dt = new TPaveText(0.5, 0.65, 0.85, 0.85, "NDC");
        const double A      = f_dt->GetParameter(0);
        const double Aerr   = f_dt->GetParError(0);
        const double tau    = f_dt->GetParameter(1);
        const double tauerr = f_dt->GetParError(1);
        const double c      = f_dt->GetParameter(2);
        const double cerr   = f_dt->GetParError(2);
        pt_dt->SetFillStyle(0);
        pt_dt->SetBorderSize(0);
        pt_dt->SetTextAlign(12);
        pt_dt->SetTextFont(42);
        pt_dt->SetTextSize(0.035);
        pt_dt->AddText("#bf{Fit:} A #upoint exp(-x / #tau) + c");
        pt_dt->AddText(Form("A = %.3g #pm %.3g", A, Aerr));
        pt_dt->AddText(Form("#tau = %.3g #pm %.3g #mus", tau, tauerr));
        pt_dt->AddText(Form("c = %.3g #pm %.3g", c, cerr));
        TCanvas* c_dt = plot_basic(h_dt, "HIST");
        f_dt->Draw("SAME");
        pt_dt->Draw("SAME");
        c_dt->Update();

        // ============================================================================================
        // Prompt-Delayed distance
        // ============================================================================================

        h_dr->SetStats(false);
        pimp_my_line(h_dr, LineConfig{.style = kSolid, .width = 3, .color = kBlue});
        pimp_my_axis(h_dr->GetXaxis(), AxisConfig{.maxdigits = 3});
        pimp_my_name(h_dr->GetXaxis(), NameConfig{.title = "d_{p-d} (mm)"});
        h_dr->GetXaxis()->CenterTitle(true);
        pimp_my_axis(h_dr->GetYaxis(), AxisConfig{.maxdigits = 3, .title = {.offset = 1.25}});
        pimp_my_name(h_dr->GetYaxis(), NameConfig{.title = "Entries"});
        h_dr->GetYaxis()->CenterTitle(true);
        TCanvas* c_dr = plot_basic(h_dr, "HIST");

        // ============================================================================================
        // Prompt vertex position
        // ============================================================================================

        h_rho_z_p->SetStats(false);
        pimp_my_axis(h_rho_z_p->GetXaxis(), AxisConfig{.maxdigits = 3});
        pimp_my_name(h_rho_z_p->GetXaxis(), NameConfig{.title = "#rho (mm)"});
        h_rho_z_p->GetXaxis()->CenterTitle(true);
        pimp_my_axis(h_rho_z_p->GetYaxis(), AxisConfig{.maxdigits = 3, .title = {.offset = 1.25}});
        pimp_my_name(h_rho_z_p->GetYaxis(), NameConfig{.title = "z (mm)"});
        h_rho_z_p->GetYaxis()->CenterTitle(true);
        TCanvas* c_rho_z_p = plot_basic(h_rho_z_p, "COLZ");

        // ============================================================================================
        // Delayed vertex position
        // ============================================================================================

        h_rho_z_d->SetStats(false);
        pimp_my_axis(h_rho_z_d->GetXaxis(), AxisConfig{.maxdigits = 3});
        pimp_my_name(h_rho_z_d->GetXaxis(), NameConfig{.title = "#rho (mm)"});
        h_rho_z_d->GetXaxis()->CenterTitle(true);
        pimp_my_axis(h_rho_z_d->GetYaxis(), AxisConfig{.maxdigits = 3, .title = {.offset = 1.25}});
        pimp_my_name(h_rho_z_d->GetYaxis(), NameConfig{.title = "z (mm)"});
        h_rho_z_d->GetYaxis()->CenterTitle(true);
        TCanvas* c_rho_z_d = plot_basic(h_rho_z_d, "COLZ");
    
    }

protected:

    std::set<ibd> m_ibds;

};

#endif // ANALYSIS_IBD_ANALYSIS_HPP_