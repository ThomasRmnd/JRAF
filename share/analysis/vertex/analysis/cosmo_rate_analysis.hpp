#ifndef ANALYSIS_COSMO_RATE_ANALYSIS_HPP_
#define ANALYSIS_COSMO_RATE_ANALYSIS_HPP_

#include <TF1.h>
#include <TFitResult.h>

#include "analysis/basic_analysis.hpp"

class cosmo_rate_analysis : public basic_analysis {

    // Regarder nombre cosmo en fonction du nombre de neutron ==> est-ce que ça suit une loi de poisson

public:

    cosmo_rate_analysis(const std::string& name, const std::string& filepath, const std::string& suffix) :
        basic_analysis{name, filepath, suffix}
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

        return (nb_multu_veto == 0ul);
    }

    bool process() override {
        ibd v{m_nav->prompt, m_nav->delayed};

        std::vector<physical_muon> phy_mu;
        std::vector<muon_data_association> muon_data;
        for (std::size_t k = 0ul; k < m_nav->method_mu.size(); ++k) {
            timestamp ts_mu{m_nav->sec_mu[k], m_nav->nsec_mu[k]};
            std::vector<physical_muon>::iterator it = std::find_if(
                phy_mu.begin(), phy_mu.end(),
                [ts_mu](const physical_muon& mu) {
                    timestamp diff = ts_mu - mu.ts;
                    return timestamp{0, -1000} < diff && diff < timestamp{0, 1000};
                }
            );
            if (it == phy_mu.end()) {
                phy_mu.push_back({ts_mu, {k}});
            }
            else {
                it->indices.push_back(k);
            }
        }

        for (const physical_muon& mu : phy_mu) {
            timestamp diff = v.prompt.ts - mu.ts;
            if (diff < timestamp{0, 5000000} || timestamp{0, 1200000000} < diff) continue;

            muon_data_association assoc;
            assoc.dt = timestamp_to_double(diff);

            int neutron_count = 0;
            for (std::size_t k = 0ul; k < m_nav->sec_n.size(); ++k) {
                double e_n = m_nav->e_n[k];
                if (e_n < 2.0 || 2.5 < e_n) continue;
                timestamp ts_n{m_nav->sec_n[k], m_nav->nsec_n[k]};
                if (ts_n < mu.ts + timestamp{0, 20000} || mu.ts + timestamp{0, 2000000} < ts_n) continue;
                ++neutron_count;
            }
            assoc.neutron_count = neutron_count;

            for (std::size_t idx : mu.indices) {
                const std::string& method = m_nav->method_mu[idx];
                if (method != "CdWpTtChi2" && method != "Tt") continue;
                vec3 pos_mu{m_nav->posx_mu[idx], m_nav->posy_mu[idx], m_nav->posz_mu[idx]};
                vec3 dir_mu{m_nav->dirx_mu[idx], m_nav->diry_mu[idx], m_nav->dirz_mu[idx]};
                double d = mag(cross(dir_mu, v.prompt.pos - pos_mu));

                if (method == "CdWpTtChi2") assoc.dlat_cdwp.push_back(d);
                else assoc.dlat_tt.push_back(d);
            }
            muon_data.push_back(std::move(assoc));
        }
        m_ibds_to_mu[v] = muon_data;

        return true;
    }

    void result() override {
        std::unordered_map<int, TH1D*> cosmo_rate_with_at_least_n_neutron;
        std::unordered_map<int, TH1D*> cosmo_rate_with_n_neutron;
        for (int k = 0; k < 10; ++k) {
            cosmo_rate_with_at_least_n_neutron[k] = new TH1D(Form("h_cosmo_rate_with_at_least_%i_neutron", k), Form("Cosmo rate with at least %i neutron", k), 120, 0.0, 1.2);
            cosmo_rate_with_n_neutron[k] = new TH1D(Form("h_cosmo_rate_with_%i_neutron", k), Form("Cosmo rate with %i neutron", k), 120, 0.0, 1.2);
        }
        for (const std::pair<ibd, std::vector<muon_data_association>>& val : m_ibds_to_mu) {
            const std::vector<muon_data_association>& muon_data = val.second;
            for (const muon_data_association& assoc : muon_data) {
                for (auto& [n, h] : cosmo_rate_with_at_least_n_neutron) {
                    if (assoc.neutron_count >= n) {
                        h->Fill(assoc.dt);
                    }
                }
                for (auto& [n, h] : cosmo_rate_with_n_neutron) {
                    if (assoc.neutron_count == n) {
                        h->Fill(assoc.dt);
                    }
                }
            }
        }

        for (auto& [n, h] : cosmo_rate_with_at_least_n_neutron) {
            fit_and_plot_cosmo_rate_with_neutron(h);
        }
        for (auto& [n, h] : cosmo_rate_with_n_neutron) {
            fit_and_plot_cosmo_rate_with_neutron(h);
        }

        TH2D* h_d_mu2p_cdwp_vs_dt_mu2p = new TH2D("h_d_mu2p_cdwp_vs_dt_mu2p", "Cosmo time vs distance", 120, 0.0, 1.5, 100, 0.0, 40000.0);
        TH2D* h_d_mu2p_tt_vs_dt_mu2p = new TH2D("h_d_mu2p_tt_vs_dt_mu2p", "Cosmo time vs distance", 120, 0.0, 1.5, 100, 0.0, 40000.0);
        for (const std::pair<ibd, std::vector<muon_data_association>>& val : m_ibds_to_mu) {
            const std::vector<muon_data_association>& muon_data = val.second;
            for (const muon_data_association& assoc : muon_data) {
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

private:

    struct physical_muon {

        timestamp ts;
        std::vector<std::size_t> indices;

    };

    struct muon_data_association {

        double dt;
        int neutron_count;
        std::vector<double> dlat_cdwp;
        std::vector<double> dlat_tt;

    };

    std::map<ibd, std::vector<muon_data_association>> m_ibds_to_mu;

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

        TFitResultPtr res = h->Fit(f, "RS");
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

        double p1 = f->GetParameter(1);
        double p2 = f->GetParameter(2);
        double sp1 = f->GetParError(1);
        double sp2 = f->GetParError(2);
    
        double covp1p2 = res->CovMatrix(1, 2);

        double exp_term = std::exp(-time_window / p2);
        double nLiHe = (p1 * p2 * (1 - exp_term)) / binning;

        double df_dp1 = (p2 * (1 - exp_term)) / binning;
        double df_dp2 = (p1 / binning) * (1 - exp_term * (1 + time_window / p2));

        double nLiHe_err = std::sqrt(std::pow(df_dp1 * sp1, 2.0) + std::pow(df_dp2 * sp2, 2.0) + 2.0 * df_dp1 * df_dp2 * covp1p2);

        std::cout << "Fit Results for " << h->GetName() << ":\n";
        std::cout << "nIBD = " << f->GetParameter(0) * time_window / binning << " +/- " << f->GetParError(0) * time_window / binning << '\n';
        std::cout << "nLiHe = " << nLiHe << " +/- " << nLiHe_err << std::endl;
    }

};

#endif // ANALYSIS_COSMO_RATE_ANALYSIS_HPP_