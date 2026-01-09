#ifndef ANALYSIS_COSMO_SHAPE_ANALYSIS_HPP_
#define ANALYSIS_COSMO_SHAPE_ANALYSIS_HPP_

#include "analysis/basic_analysis.hpp"

class cosmo_shape_analysis : public basic_analysis {

    // Essayer de voir si on peut, par le calcul, retrouver le nombre théorique de cosmo en prenant en compte les inéfficacité
    // Essayer de faire la même chose mais avec les autres méthodes de reconstruction, i.e. CdClassify et WpClassify

public:

    cosmo_shape_analysis(
        const std::string& name, 
        const std::string& filepath, const std::string& suffix, 
        const std::string& recname,
        const timestamp& sig_low, const timestamp& sig_high, 
        const timestamp& bkg_low, const timestamp& bkg_high, 
        double radius
    ) :
        basic_analysis{name, filepath, suffix},
        m_recname{recname},
        m_ts_sig_low{sig_low},
        m_ts_sig_high{sig_high},
        m_ts_bkg_low{bkg_low},
        m_ts_bkg_high{bkg_high},
        m_radius{radius}
    {}

    bool selection() override {
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

        m_dlat_mu2p.clear();
        m_dlat_mu2d.clear();
        m_dt_mu2p.clear();
        m_dt_mu2d.clear();
        m_is_sig.clear();

        for (std::size_t k = 0ul; k < m_nav->method_mu.size(); ++k) {
            if (m_nav->method_mu[k] != m_recname) continue;
            timestamp ts_mu{m_nav->sec_mu[k], m_nav->nsec_mu[k]};
            bool is_in_bkg = (
                ts_mu + m_ts_bkg_low < m_nav->prompt.ts && m_nav->prompt.ts < ts_mu + m_ts_bkg_high &&
                ts_mu + m_ts_bkg_low < m_nav->delayed.ts && m_nav->delayed.ts < ts_mu + m_ts_bkg_high
            );
            bool is_in_sig = (
                ts_mu + m_ts_sig_low < m_nav->prompt.ts && m_nav->prompt.ts < ts_mu + m_ts_sig_high &&
                ts_mu + m_ts_sig_low < m_nav->delayed.ts && m_nav->delayed.ts < ts_mu + m_ts_sig_high
            );
            if (!is_in_bkg && !is_in_sig) continue;
            vec3 pos_mu{m_nav->posx_mu[k], m_nav->posy_mu[k], m_nav->posz_mu[k]};
            vec3 dir_mu{m_nav->dirx_mu[k], m_nav->diry_mu[k], m_nav->dirz_mu[k]};
            double d_mu2p = mag(cross(dir_mu, m_nav->prompt.pos - pos_mu));
            double d_mu2d = mag(cross(dir_mu, m_nav->delayed.pos - pos_mu));
            if (m_radius < d_mu2p || m_radius < d_mu2d) continue;
            m_dlat_mu2p.push_back(d_mu2p);
            m_dlat_mu2d.push_back(d_mu2d);
            m_dt_mu2p.push_back(timestamp_to_double(m_nav->prompt.ts - ts_mu));
            m_dt_mu2d.push_back(timestamp_to_double(m_nav->delayed.ts - ts_mu));
            m_is_sig.push_back(is_in_sig);
        }

        return (nb_multu_veto == 0ul && !m_is_sig.empty());
    }

    bool process() override {
        for (std::size_t k = 0ul; k < m_is_sig.size(); ++k) {
            if (m_is_sig[k]) {
                m_cosmos_sig.insert({m_nav->prompt, m_nav->delayed, m_dlat_mu2p[k], m_dlat_mu2d[k], m_dt_mu2p[k], m_dt_mu2d[k]});
            }
            else {
                m_cosmos_bkg.insert({m_nav->prompt, m_nav->delayed, m_dlat_mu2p[k], m_dlat_mu2d[k], m_dt_mu2p[k], m_dt_mu2d[k]});
            }
        }
        return true;
    }

    void result() override {
        std::vector<cosmogenic> cosmos_bkg;
        cosmos_bkg.reserve(m_cosmos_bkg.size());
        for (std::set<cosmogenic>::const_iterator it = m_cosmos_bkg.begin(); it != m_cosmos_bkg.end(); ++it) {
            cosmos_bkg.push_back(*it);
        }
        std::vector<cosmogenic> cosmos_sig;
        cosmos_sig.reserve(m_cosmos_sig.size());
        for (std::set<cosmogenic>::const_iterator it = m_cosmos_sig.begin(); it != m_cosmos_sig.end(); ++it) {
            cosmos_sig.push_back(*it);
        }

        TH1D* h_e_p_cosmo_bkg = make_prompt_energy_plot(Form("h_e_p_cosmo_bkg__%s", m_name.c_str()), Form("Prompt energy (Cosmo bkg) {%s}", m_name.c_str()), cosmos_bkg);
        TH1D* h_e_p_cosmo_sig = make_prompt_energy_plot(Form("h_e_p_cosmo_sig__%s", m_name.c_str()), Form("Prompt energy (Cosmo sig) {%s}", m_name.c_str()), cosmos_sig);
        TH1D* h_e_p_cosmo_diff = make_prompt_energy_plot(Form("h_e_p_cosmo_diff__%s", m_name.c_str()), Form("Prompt energy (Cosmo sig - Cosmo bkg) {%s}", m_name.c_str()), std::vector<cosmogenic>{});

        TH1D* h_e_d_cosmo_bkg = make_delayed_energy_plot(Form("h_e_d_cosmo_bkg__%s", m_name.c_str()), Form("Delayed energy (Cosmo bkg) {%s}", m_name.c_str()), cosmos_bkg);
        TH1D* h_e_d_cosmo_sig = make_delayed_energy_plot(Form("h_e_d_cosmo_sig__%s", m_name.c_str()), Form("Delayed energy (Cosmo sig) {%s}", m_name.c_str()), cosmos_sig);

        TH1D* h_dt_cosmo_bkg = make_prompt_delayed_time_plot(Form("h_dt_cosmo_bkg__%s", m_name.c_str()), Form("Prompt-Delayed time difference (Cosmo bkg) {%s}", m_name.c_str()), cosmos_bkg); 
        TH1D* h_dt_cosmo_sig = make_prompt_delayed_time_plot(Form("h_dt_cosmo_sig__%s", m_name.c_str()), Form("Prompt-Delayed time difference (Cosmo sig) {%s}", m_name.c_str()), cosmos_sig);

        TH1D* h_dr_cosmo_bkg = make_prompt_delayed_distance_plot(Form("h_dr_cosmo_bkg__%s", m_name.c_str()), Form("Prompt-Delayed distance (Cosmo bkg) {%s}", m_name.c_str()), cosmos_bkg);
        TH1D* h_dr_cosmo_sig = make_prompt_delayed_distance_plot(Form("h_dr_cosmo_sig__%s", m_name.c_str()), Form("Prompt-Delayed distance (Cosmo sig) {%s}", m_name.c_str()), cosmos_sig);

        TH2D* h_rho_z_p_cosmo_bkg = make_prompt_spatial_plot(Form("h_rho_z_p_cosmo_bkg__%s", m_name.c_str()), Form("Prompt vertex distribution (Cosmo bkg) {%s}", m_name.c_str()), cosmos_bkg);
        TH2D* h_rho_z_d_cosmo_bkg = make_delayed_spatial_plot(Form("h_rho_z_d_cosmo_bkg__%s", m_name.c_str()), Form("Delayed vertex distribution (Cosmo bkg) {%s}", m_name.c_str()), cosmos_bkg);
        TH2D* h_rho_z_p_cosmo_sig = make_prompt_spatial_plot(Form("h_rho_z_p_cosmo_sig__%s", m_name.c_str()), Form("Prompt vertex distribution (Cosmo sig) {%s}", m_name.c_str()), cosmos_sig);
        TH2D* h_rho_z_d_cosmo_sig = make_delayed_spatial_plot(Form("h_rho_z_d_cosmo_sig__%s", m_name.c_str()), Form("Delayed vertex distribution (Cosmo sig) {%s}", m_name.c_str()), cosmos_sig);

        h_e_p_cosmo_diff->Add(h_e_p_cosmo_sig, h_e_p_cosmo_bkg, 1.0, -1.0);

        // ============================================================================================
        // Cosmo bkg - Prompt energy
        // ============================================================================================

        pimp_my_histogram(h_e_p_cosmo_bkg, kSolid, 3, kBlue, 1.0);
        plot_basic(h_e_p_cosmo_bkg);

        // ============================================================================================
        // Cosmo bkg - Delayed energy
        // ============================================================================================

        TCanvas* c_e_d_cosmo_bkg = new TCanvas("c_e_d_cosmo_bkg", "Delayed energy (Cosmo bkg)", 1000, 1000);
        c_e_d_cosmo_bkg->cd();

        pimp_my_histogram(h_e_d_cosmo_bkg, kSolid, 3, kBlue, 1.0);
        h_e_d_cosmo_bkg->Draw();

        c_e_d_cosmo_bkg->Update();

        // ============================================================================================
        // Cosmo bkg - Prompt-Delayed time difference
        // ============================================================================================

        TCanvas* c_dt_cosmo_bkg = new TCanvas("c_dt_cosmo_bkg", "Prompt-Delayed time difference (Cosmo bkg)", 1000, 1000);
        c_dt_cosmo_bkg->cd();

        pimp_my_histogram(h_dt_cosmo_bkg, kSolid, 3, kBlue, 1.0);
        h_dt_cosmo_bkg->Draw();

        c_dt_cosmo_bkg->Update();

        // ============================================================================================
        // Cosmo bkg - Prompt-Delayed distance
        // ============================================================================================

        TCanvas* c_dr_cosmo_bkg = new TCanvas("c_dr_cosmo_bkg", "Prompt-Delayed distance (Cosmo bkg)", 1000, 1000);
        c_dr_cosmo_bkg->cd();

        h_dr_cosmo_bkg->Draw();

        c_dr_cosmo_bkg->Update();

        // ============================================================================================
        // Cosmo bkg - Prompt vertex position
        // ============================================================================================

        TCanvas* c_rho_z_p_cosmo_bkg = new TCanvas("c_rho_z_p_cosmo_bkg", "Prompt vertex distribution (Cosmo bkg)", 1000, 1000);
        c_rho_z_p_cosmo_bkg->cd();

        h_rho_z_p_cosmo_bkg->Draw();

        c_rho_z_p_cosmo_bkg->Update();

        // ============================================================================================
        // Cosmo bkg - Delayed vertex position
        // ============================================================================================

        TCanvas* c_rho_z_d_cosmo_bkg = new TCanvas("c_rho_z_d_cosmo_bkg", "Delayed vertex distribution (Cosmo bkg)", 1000, 1000);
        c_rho_z_d_cosmo_bkg->cd();

        h_rho_z_d_cosmo_bkg->Draw();

        c_rho_z_d_cosmo_bkg->Update();

        // ============================================================================================
        // Cosmo sig - Prompt energy
        // ============================================================================================

        TCanvas* c_e_p_cosmo_sig = new TCanvas("c_e_p_cosmo_sig", "Prompt energy (Cosmo sig)", 1000, 1000);
        c_e_p_cosmo_sig->cd();

        h_e_p_cosmo_sig->SetLineWidth(3);
        h_e_p_cosmo_sig->SetLineStyle(kSolid);
        h_e_p_cosmo_sig->SetLineColorAlpha(kBlue, 1.0);

        h_e_p_cosmo_sig->Draw();

        c_e_p_cosmo_sig->Update();

        // ============================================================================================
        // Cosmo sig - Delayed energy
        // ============================================================================================

        TCanvas* c_e_d_cosmo_sig = new TCanvas("c_e_d_cosmo_sig", "Delayed energy (Cosmo sig)", 1000, 1000);
        c_e_d_cosmo_sig->cd();
    
        h_e_d_cosmo_sig->SetLineWidth(3);
        h_e_d_cosmo_sig->SetLineStyle(kSolid);
        h_e_d_cosmo_sig->SetLineColorAlpha(kBlue, 1.0);

        h_e_d_cosmo_sig->Draw();

        c_e_d_cosmo_sig->Update();

        // ============================================================================================
        // Cosmo sig - Prompt-Delayed time difference
        // ============================================================================================

        TCanvas* c_dt_cosmo_sig = new TCanvas("c_dt_cosmo_sig", "Prompt-Delayed time difference (Cosmo sig)", 1000, 1000);
        c_dt_cosmo_sig->cd();

        h_dt_cosmo_sig->SetLineWidth(3);
        h_dt_cosmo_sig->SetLineStyle(kSolid);
        h_dt_cosmo_sig->SetLineColorAlpha(kBlue, 1.0);

        h_dt_cosmo_sig->Draw();

        c_dt_cosmo_sig->Update();

        // ============================================================================================
        // Cosmo sig - Prompt-Delayed distance
        // ============================================================================================

        TCanvas* c_dr_cosmo_sig = new TCanvas("c_dr_cosmo_sig", "Prompt-Delayed distance (Cosmo sig)", 1000, 1000);
        c_dr_cosmo_sig->cd();

        h_dr_cosmo_sig->Draw();

        c_dr_cosmo_sig->Update();

        // ============================================================================================
        // Cosmo sig - Prompt vertex position
        // ============================================================================================

        TCanvas* c_rho_z_p_cosmo_sig = new TCanvas("c_rho_z_p_cosmo_sig", "Prompt vertex distribution (Cosmo sig)", 1000, 1000);
        c_rho_z_p_cosmo_sig->cd();

        h_rho_z_p_cosmo_sig->Draw();

        c_rho_z_p_cosmo_sig->Update();

        // ============================================================================================
        // Cosmo sig - Delayed vertex position
        // ============================================================================================

        TCanvas* c_rho_z_d_cosmo_sig = new TCanvas("c_rho_z_d_cosmo_sig", "Delayed vertex distribution (Cosmo sig)", 1000, 1000);
        c_rho_z_d_cosmo_sig->cd();

        h_rho_z_d_cosmo_sig->Draw();

        c_rho_z_d_cosmo_sig->Update();

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

private:

    std::string m_recname;
    timestamp m_ts_sig_low;
    timestamp m_ts_sig_high;
    timestamp m_ts_bkg_low;
    timestamp m_ts_bkg_high;
    double m_radius;

    std::vector<double> m_dlat_mu2p;
    std::vector<double> m_dlat_mu2d;
    std::vector<double> m_dt_mu2p;
    std::vector<double> m_dt_mu2d;
    std::vector<bool> m_is_sig;

    std::set<cosmogenic> m_cosmos_bkg, m_cosmos_sig;

};

#endif // ANALYSIS_COSMO_SHAPE_ANALYSIS_HPP_