#ifndef ANALYSIS_COSMO_SHAPE_NEUTRON_ANALYSIS_HPP_
#define ANALYSIS_COSMO_SHAPE_NEUTRON_ANALYSIS_HPP_

#include <set>

#include <TFile.h>

#include "analysis/basic_analysis.hpp"
#include "utils/plot.hpp"

class cosmo_shape_neutron_analysis : public basic_analysis {

public:

    cosmo_shape_neutron_analysis(
        const std::string& name, 
        const std::string& filepath, const std::string& suffix,
        const timestamp& sig_low, const timestamp& sig_high, 
        const timestamp& bkg_low, const timestamp& bkg_high, 
        double radius
    ) :
        basic_analysis{name, filepath, suffix},
        m_ts_sig_low{sig_low},
        m_ts_sig_high{sig_high},
        m_ts_bkg_low{bkg_low},
        m_ts_bkg_high{bkg_high},
        m_radius{radius}
    {}

    virtual ~cosmo_shape_neutron_analysis() override = default;

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

        m_d_neu2p.clear();
        m_d_neu2d.clear();
        m_dt_neu2p.clear();
        m_dt_neu2d.clear();
        m_is_sig.clear();

        for (std::size_t k = 0ul; k < m_nav->e_n.size(); ++k) {
            if (m_nav->e_n[k] < 1.5 || 20.0 < m_nav->e_n[k]) continue;
            timestamp ts_n{m_nav->sec_n[k], m_nav->nsec_n[k]};
            bool is_in_bkg = (
                ts_n + m_ts_bkg_low < m_nav->prompt.ts && m_nav->prompt.ts < ts_n + m_ts_bkg_high &&
                ts_n + m_ts_bkg_low < m_nav->delayed.ts && m_nav->delayed.ts < ts_n + m_ts_bkg_high
            );
            bool is_in_sig = (
                ts_n + m_ts_sig_low < m_nav->prompt.ts && m_nav->prompt.ts < ts_n + m_ts_sig_high &&
                ts_n + m_ts_sig_low < m_nav->delayed.ts && m_nav->delayed.ts < ts_n + m_ts_sig_high
            );
            if (!is_in_bkg && !is_in_sig) continue;
            vec3 pos_n{m_nav->posx_n[k], m_nav->posy_n[k], m_nav->posz_n[k]};
            double d_neu2p = mag(m_nav->prompt.pos - pos_n);
            double d_neu2d = mag(m_nav->delayed.pos - pos_n);
            if (m_radius < d_neu2p || m_radius < d_neu2d) continue;
            m_d_neu2p.push_back(d_neu2p);
            m_d_neu2d.push_back(d_neu2d);
            m_dt_neu2p.push_back(timestamp_to_double(m_nav->prompt.ts - ts_n));
            m_dt_neu2d.push_back(timestamp_to_double(m_nav->delayed.ts - ts_n));
            m_is_sig.push_back(is_in_sig);
        }

        return (nb_multu_veto == 0ul && !m_is_sig.empty());
    }

    bool process() override {
        for (std::size_t k = 0ul; k < m_is_sig.size(); ++k) {
            // cosmo has members dlat_mu2p, dt_mu2p, etc. but we will use them as d_neu2p, dt_neu2p, etc. 
            if (m_is_sig[k]) {
                m_cosmos_sig.insert({m_nav->run_id, m_nav->prompt, m_nav->delayed, m_d_neu2p[k], m_d_neu2d[k], m_dt_neu2p[k], m_dt_neu2d[k]});
            }
            else {
                m_cosmos_bkg.insert({m_nav->run_id, m_nav->prompt, m_nav->delayed, m_d_neu2p[k], m_d_neu2d[k], m_dt_neu2p[k], m_dt_neu2d[k]});
            }
        }
        return true;
    }

    bool save() override {
        TFile* f = TFile::Open(Form("%s.root", m_name.c_str()), "RECREATE");
        if (!f) {
            std::cerr << "Cannot open file " << m_name << ".root for writing\n";
            return false;
        }
        TTree* t_bkg = new TTree("background_events", "Cosmogenic depleted region");
        TTree* t_sig = new TTree("signal_events", "Cosmogenic enriched region");
        if (!t_bkg || !t_sig) {
            std::cerr << "Cannot create tree background or signal\n";
            return false;
        }
        int run_id;
        vec3 pos_p, pos_d;
        timestamp ts_p, ts_d;
        double e_p, e_d;
        double dlat_mu2p, dlat_mu2d; 
        double dt_mu2p, dt_mu2d;
        
        t_bkg->Branch("run_id", &run_id);
        t_bkg->Branch("posx_p", &pos_p.x);
        t_bkg->Branch("posy_p", &pos_p.y);
        t_bkg->Branch("posz_p", &pos_p.z);
        t_bkg->Branch("sec_p", &ts_p.sec);
        t_bkg->Branch("nsec_p", &ts_p.nsec);
        t_bkg->Branch("e_p", &e_p);
        t_bkg->Branch("dlat_mu2p", &dlat_mu2p);
        t_bkg->Branch("dt_mu2p", &dt_mu2p);
        t_bkg->Branch("posx_d", &pos_d.x);
        t_bkg->Branch("posy_d", &pos_d.y);
        t_bkg->Branch("posz_d", &pos_d.z);
        t_bkg->Branch("sec_d", &ts_d.sec);
        t_bkg->Branch("nsec_d", &ts_d.nsec);
        t_bkg->Branch("e_d", &e_d);
        t_bkg->Branch("dlat_mu2d", &dlat_mu2d);
        t_bkg->Branch("dt_mu2d", &dt_mu2d);

        for (std::set<cosmogenic>::const_iterator it = m_cosmos_bkg.begin(); it != m_cosmos_bkg.end(); ++it) {
            run_id = it->run_id;
            pos_p = it->prompt.pos;
            pos_d = it->delayed.pos;
            ts_p = it->prompt.ts;
            ts_d = it->delayed.ts;
            e_p = it->prompt.e;
            e_d = it->delayed.e;
            dlat_mu2p = it->dlat_mu2p;
            dlat_mu2d = it->dlat_mu2d;
            dt_mu2p = it->dt_mu2p;
            dt_mu2d = it->dt_mu2d;
            t_bkg->Fill();
        }

        t_sig->Branch("run_id", &run_id);
        t_sig->Branch("posx_p", &pos_p.x);
        t_sig->Branch("posy_p", &pos_p.y);
        t_sig->Branch("posz_p", &pos_p.z);
        t_sig->Branch("sec_p", &ts_p.sec);
        t_sig->Branch("nsec_p", &ts_p.nsec);
        t_sig->Branch("e_p", &e_p);
        t_sig->Branch("dlat_mu2p", &dlat_mu2p);
        t_sig->Branch("dt_mu2p", &dt_mu2p);
        t_sig->Branch("posx_d", &pos_d.x);
        t_sig->Branch("posy_d", &pos_d.y);
        t_sig->Branch("posz_d", &pos_d.z);
        t_sig->Branch("sec_d", &ts_d.sec);
        t_sig->Branch("nsec_d", &ts_d.nsec);
        t_sig->Branch("e_d", &e_d);
        t_sig->Branch("dlat_mu2d", &dlat_mu2d);
        t_sig->Branch("dt_mu2d", &dt_mu2d);

        for (std::set<cosmogenic>::const_iterator it = m_cosmos_sig.begin(); it != m_cosmos_sig.end(); ++it) {
            run_id = it->run_id;
            pos_p = it->prompt.pos;
            pos_d = it->delayed.pos;
            ts_p = it->prompt.ts;
            ts_d = it->delayed.ts;
            e_p = it->prompt.e;
            e_d = it->delayed.e;
            dlat_mu2p = it->dlat_mu2p;
            dlat_mu2d = it->dlat_mu2d;
            dt_mu2p = it->dt_mu2p;
            dt_mu2d = it->dt_mu2d;
            t_sig->Fill();
        }

        f->cd();
        t_bkg->Write();
        t_sig->Write();
        f->Close();
        return true;
    }

    void result() override {
        std::vector<cosmogenic> cosmos_bkg(m_cosmos_bkg.begin(), m_cosmos_bkg.end());
        std::vector<cosmogenic> cosmos_sig(m_cosmos_sig.begin(), m_cosmos_sig.end());

        TH1D* h_e_p_cosmo_bkg = make_normal_prompt_energy_plot(Form("h_e_p_cosmo_bkg__%s", m_name.c_str()), Form("Prompt energy (Cosmo bkg) {%s}", m_name.c_str()), cosmos_bkg);
        TH1D* h_e_p_cosmo_sig = make_normal_prompt_energy_plot(Form("h_e_p_cosmo_sig__%s", m_name.c_str()), Form("Prompt energy (Cosmo sig) {%s}", m_name.c_str()), cosmos_sig);
        TH1D* h_e_p_cosmo_diff = make_normal_prompt_energy_plot(Form("h_e_p_cosmo_diff__%s", m_name.c_str()), Form("Prompt energy (Cosmo sig - Cosmo bkg) {%s}", m_name.c_str()), std::vector<cosmogenic>{});
        h_e_p_cosmo_diff->Add(h_e_p_cosmo_sig, h_e_p_cosmo_bkg, 1.0, -1.0);

        TH1D* h_e_d_cosmo_bkg = make_delayed_energy_plot(Form("h_e_d_cosmo_bkg__%s", m_name.c_str()), Form("Delayed energy (Cosmo bkg) {%s}", m_name.c_str()), cosmos_bkg);
        TH1D* h_e_d_cosmo_sig = make_delayed_energy_plot(Form("h_e_d_cosmo_sig__%s", m_name.c_str()), Form("Delayed energy (Cosmo sig) {%s}", m_name.c_str()), cosmos_sig);

        TH1D* h_dt_cosmo_bkg = make_prompt_delayed_time_plot(Form("h_dt_cosmo_bkg__%s", m_name.c_str()), Form("Prompt-Delayed time difference (Cosmo bkg) {%s}", m_name.c_str()), cosmos_bkg); 
        TH1D* h_dt_cosmo_sig = make_prompt_delayed_time_plot(Form("h_dt_cosmo_sig__%s", m_name.c_str()), Form("Prompt-Delayed time difference (Cosmo sig) {%s}", m_name.c_str()), cosmos_sig);

        TH1D* h_dr_cosmo_bkg = make_prompt_delayed_distance_plot(Form("h_dr_cosmo_bkg__%s", m_name.c_str()), Form("Prompt-Delayed distance (Cosmo bkg) {%s}", m_name.c_str()), cosmos_bkg);
        TH1D* h_dr_cosmo_sig = make_prompt_delayed_distance_plot(Form("h_dr_cosmo_sig__%s", m_name.c_str()), Form("Prompt-Delayed distance (Cosmo sig) {%s}", m_name.c_str()), cosmos_sig);

        // TH2D* h_rho_z_p_cosmo_bkg = make_prompt_spatial_plot(Form("h_rho_z_p_cosmo_bkg__%s", m_name.c_str()), Form("Prompt vertex distribution (Cosmo bkg) {%s}", m_name.c_str()), cosmos_bkg);
        // TH2D* h_rho_z_d_cosmo_bkg = make_delayed_spatial_plot(Form("h_rho_z_d_cosmo_bkg__%s", m_name.c_str()), Form("Delayed vertex distribution (Cosmo bkg) {%s}", m_name.c_str()), cosmos_bkg);
        // TH2D* h_rho_z_p_cosmo_sig = make_prompt_spatial_plot(Form("h_rho_z_p_cosmo_sig__%s", m_name.c_str()), Form("Prompt vertex distribution (Cosmo sig) {%s}", m_name.c_str()), cosmos_sig);
        // TH2D* h_rho_z_d_cosmo_sig = make_delayed_spatial_plot(Form("h_rho_z_d_cosmo_sig__%s", m_name.c_str()), Form("Delayed vertex distribution (Cosmo sig) {%s}", m_name.c_str()), cosmos_sig);

        // Prompt energy
        TLegend* leg_e_p = new TLegend(0.58, 0.75, 0.88, 0.88);
        pimp_my_line(h_e_p_cosmo_bkg, LineConfig{.style = kSolid, .width = 3, .color = kRed});
        h_e_p_cosmo_bkg->SetTitle("Cosmo depleted region");
        pimp_my_line(h_e_p_cosmo_sig, LineConfig{.style = kSolid, .width = 3, .color = kBlue});
        h_e_p_cosmo_sig->SetTitle("Cosmo enriched region");
        pimp_my_line(h_e_p_cosmo_diff, LineConfig{.style = kSolid, .width = 3, .color = kBlack});
        h_e_p_cosmo_diff->SetTitle("Substraction");
        plot_multiple(Form("c_e_p_cosmo__%s", m_name.c_str()), "Prompt energy", {h_e_p_cosmo_bkg, h_e_p_cosmo_sig, h_e_p_cosmo_diff}, leg_e_p, "HIST");

        // Cosmo bkg - Delayed energy
        TLegend* leg_e_d = new TLegend(0.58, 0.75, 0.88, 0.88);
        pimp_my_line(h_e_d_cosmo_bkg, LineConfig{.style = kSolid, .width = 3, .color = kRed});
        h_e_d_cosmo_bkg->SetTitle("Cosmo depleted region");
        pimp_my_line(h_e_d_cosmo_sig, LineConfig{.style = kSolid, .width = 3, .color = kBlue});
        h_e_d_cosmo_sig->SetTitle("Cosmo enriched region");
        plot_multiple(Form("c_e_d_cosmo__%s", m_name.c_str()), "Delayed energy", {h_e_d_cosmo_bkg, h_e_d_cosmo_sig}, leg_e_d, "HIST");

        // Cosmo bkg - Prompt-Delayed time difference
        TLegend* leg_dt = new TLegend(0.58, 0.75, 0.88, 0.88);
        pimp_my_line(h_dt_cosmo_bkg, LineConfig{.style = kSolid, .width = 3, .color = kRed});
        h_dt_cosmo_bkg->SetTitle("Cosmo depleted region");
        pimp_my_line(h_dt_cosmo_sig, LineConfig{.style = kSolid, .width = 3, .color = kBlue});
        h_dt_cosmo_sig->SetTitle("Cosmo enriched region");
        plot_multiple(Form("c_dt_cosmo__%s", m_name.c_str()), "Prompt-Delayed time difference", {h_dt_cosmo_bkg, h_dt_cosmo_sig}, leg_dt, "HIST");

        // Cosmo bkg - Prompt-Delayed distance
        TLegend* leg_dr = new TLegend(0.58, 0.75, 0.88, 0.88);
        pimp_my_line(h_dr_cosmo_bkg, LineConfig{.style = kSolid, .width = 3, .color = kRed});
        h_dr_cosmo_bkg->SetTitle("Cosmo depleted region");
        pimp_my_line(h_dr_cosmo_sig, LineConfig{.style = kSolid, .width = 3, .color = kBlue});
        h_dr_cosmo_sig->SetTitle("Cosmo enriched region");
        plot_multiple(Form("c_dr_cosmo__%s", m_name.c_str()), "Prompt-Delayed distance", {h_dr_cosmo_bkg, h_dr_cosmo_sig}, leg_dr, "HIST");

        // Cosmo bkg - Prompt vertex position
        // plot_basic(h_rho_z_p_cosmo_bkg, "COLZ");

        // Cosmo bkg - Delayed vertex position
        // plot_basic(h_rho_z_d_cosmo_bkg, "COLZ");

        // Cosmo sig - Prompt vertex position
        // plot_basic(h_rho_z_p_cosmo_sig, "COLZ");

        // Cosmo sig - Delayed vertex position
        // plot_basic(h_rho_z_d_cosmo_sig, "COLZ");
    }


private:

    timestamp m_ts_sig_low;
    timestamp m_ts_sig_high;
    timestamp m_ts_bkg_low;
    timestamp m_ts_bkg_high;
    double m_radius;

    std::vector<double> m_d_neu2p;
    std::vector<double> m_d_neu2d;
    std::vector<double> m_dt_neu2p;
    std::vector<double> m_dt_neu2d;
    std::vector<bool> m_is_sig;

    std::set<cosmogenic> m_cosmos_bkg, m_cosmos_sig;

};

#endif // ANALYSIS_COSMO_SHAPE_NEUTRON_ANALYSIS_HPP_