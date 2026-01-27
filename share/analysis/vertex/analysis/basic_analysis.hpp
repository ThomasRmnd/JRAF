#ifndef ANALYSIS_BASIC_ANALYSIS_HPP_
#define ANALYSIS_BASIC_ANALYSIS_HPP_

#include <set>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include "analysis/analysis.hpp"
#include "utils/plot.hpp"

class basic_analysis : public analysis_base {

public:

    basic_analysis(const std::string& name, const std::string& filepath, const std::string& suffix) :
        analysis_base{name}
    {
        std::string treename = "IBDAnalysis" + suffix;
        m_nav = navigator_manager::retrieve<basic_navigator>(filepath, treename);
        if (!m_nav->is_valid()) {
            std::cerr << "Cannot retrieve navigator of filepath " << filepath << " and treename " << treename << '\n';
            return;
        }
    }

    ~basic_analysis() override = default;

    virtual std::shared_ptr<navigator_base> navigator() const override { return m_nav; }

    virtual bool selection() override {
        if (m_nav->meta_prompt.stdt > 200.0 || m_nav->meta_delayed.stdt > 200.0) return false; // Flasher cut
        if (mag(m_nav->prompt.pos) > 16500.0) return false; // Fiducial cut
        if ((m_nav->prompt.pos.z < -15500.0 || 15500 < m_nav->prompt.pos.z) && std::sqrt(m_nav->prompt.pos.x * m_nav->prompt.pos.x + m_nav->prompt.pos.y * m_nav->prompt.pos.y) < 3000.0) return false; // Chimney cut
        if (m_nav->prompt.e < 0.7 || 12.0 < m_nav->prompt.e) return false; // Prompt energy cut
        if (m_nav->delayed.e < 2.0 || 2.5 < m_nav->delayed.e) return false; // Delayed energy cut

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < m_nav->e_mult.size(); ++k) {
            if (m_nav->e_mult[k] < 2.0 || 12.0 < m_nav->e_mult[k]) continue;
            timestamp ts_mult{m_nav->sec_mult[k], m_nav->nsec_mult[k]};
            vec3 pos_mult{m_nav->posx_mult[k], m_nav->posy_mult[k], m_nav->posz_mult[k]};
            // if (pos_mult.Mag() > 17700.0) continue;
            // if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_mult < m_nav->prompt.ts - timestamp{0, 1000000} || m_nav->delayed.ts + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }

        std::size_t nb_neutron_veto = 0ul;
        for (std::size_t k = 0ul; k < m_nav->e_n.size(); ++k) {
            if (m_nav->e_n[k] < 1.5 || 20.0 < m_nav->e_n[k]) continue;
            timestamp ts_n{m_nav->sec_n[k], m_nav->nsec_n[k]};
            vec3 pos_n{m_nav->posx_n[k], m_nav->posy_n[k], m_nav->posz_n[k]};
            // if (pos_n.Mag() > 17700.0) continue;
            if (mag(m_nav->prompt.pos - pos_n) > 4000.0 || mag(m_nav->delayed.pos - pos_n) > 4000.0) continue;
            if (m_nav->prompt.ts < ts_n + timestamp{0, 20000} || ts_n + timestamp{0, 1200000000} < m_nav->prompt.ts) continue;
            if (m_nav->delayed.ts < ts_n + timestamp{0, 20000} || ts_n + timestamp{0, 1200000000} < m_nav->delayed.ts) continue;
            ++nb_neutron_veto;
        }

        return (nb_neutron_veto == 0ul && nb_multu_veto == 0ul);
    }

    virtual bool process() override {
        m_ibds.insert({m_nav->prompt, m_nav->delayed});
        return true;
    }

    virtual void result() override {
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
        pimp_my_fill(h_e_p, FillConfig{.style = kSolid, .color = kAzure - 4});
        pimp_my_line(h_e_p, LineConfig{.style = kSolid, .width = 3, .color = kAzure - 2});
        pimp_my_axis(h_e_p->GetXaxis(), AxisConfig{.ndivisions = 406, .maxdigits = 3});
        pimp_my_name(h_e_p->GetXaxis(), NameConfig{.title = "E_{p} (MeV)"});
        h_e_p->GetXaxis()->CenterTitle(true);
        pimp_my_axis(h_e_p->GetYaxis(), AxisConfig{.maxdigits = 3, .title = {.offset = 1.25}});
        pimp_my_name(h_e_p->GetYaxis(), NameConfig{.title = "Entries"});
        h_e_p->GetYaxis()->CenterTitle(true);
        TCanvas* c_e_p = plot_basic(h_e_p, "E2");

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
        TCanvas* c_e_d = plot_basic(h_e_d, "HIST");

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
        // Write the fit result in a text box with written:
        // Fit: A * exp(-x / \tau) + c
        // A = [...] +/- [...]
        // \tau = [...] +/- [...] us
        // c = [...] +/- [...]
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

    std::shared_ptr<basic_navigator> m_nav;

    std::set<ibd> m_ibds;

};

using ibd_analysis = basic_analysis;

#endif // ANALYSIS_BASIC_ANALYSIS_HPP_