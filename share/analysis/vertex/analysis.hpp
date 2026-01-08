#ifndef ANALYSIS_HPP_
#define ANALYSIS_HPP_

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>

#include "utils/chain_reader.hpp"
#include "utils/event.hpp"
#include "utils/navigator.hpp"
#include "utils/numpy.hpp"

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

class analysis_base {

public:

    virtual ~analysis_base() = default;

    virtual std::shared_ptr<navigator_base> navigator() const = 0;

    virtual bool selection() = 0;
    virtual bool process() = 0;
    
    virtual void result() = 0;

protected:

};

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

class basic_analysis : public analysis_base {

public:

    basic_analysis(const std::string& filepath, const std::string& suffix) {
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

        TFile* vanessa_file = TFile::Open("/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/ibd/summary/ReProd25C/IBD_all_reprodC.root", "READ");
        TTree* vanessa_tree = vanessa_file->Get<TTree>("events");
        std::set<VanessaIBD> vanessa_ibds_ordered;
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
            VanessaIBD vanessa_ibd;
            vanessa_ibd.run_id = run_id;
            vanessa_ibd.ts_p = timestamp{sec_p, nsec_p};
            vanessa_ibd.e_p = e_p;
            vanessa_ibd.totq_p = totq_p;
            vanessa_ibd.ts_d = timestamp{sec_d, nsec_d};
            vanessa_ibd.e_d = e_d;
            vanessa_ibd.totq_d = totq_d;
            vanessa_ibds_ordered.insert(vanessa_ibd);
        }
        std::vector<VanessaIBD> vanessa_ibds;
        vanessa_ibds.reserve(vanessa_ibds_ordered.size());
        for (std::set<VanessaIBD>::const_iterator it = vanessa_ibds_ordered.begin(); it != vanessa_ibds_ordered.end(); ++it) {
            vanessa_ibds.push_back(*it);
        }
        
        // double e_p_min = 0.7;
        // double e_p_max = 12.0;
        // double e_p_width = 0.20;
        // int e_p_nbin = std::round((e_p_max - e_p_min) / e_p_width) + 1;
        // std::vector<double> e_p_bins = np::linspace(e_p_min, e_p_max, e_p_nbin);
        std::vector<double> e_p_bins = create_custom_e_p_bins();
        TH1D* h_e_p = new TH1D("h_e_p", "Prompt energy", e_p_bins.size() - 1, e_p_bins.data());
        TH1D* h_e_p_vanessa = new TH1D("h_e_p_vanessa", "Prompt energy (Vanessa)", e_p_bins.size() - 1, e_p_bins.data());

        double e_d_min = 2.0;
        double e_d_max = 2.5;
        double e_d_width = 0.02;
        int e_d_nbin = std::round((e_d_max - e_d_min) / e_d_width) + 1;
        std::vector<double> e_d_bins = np::linspace(e_d_min, e_d_max, e_d_nbin);
        TH1D* h_e_d = new TH1D("h_e_d", "Delayed energy", e_d_bins.size() - 1, e_d_bins.data());
        TH1D* h_e_d_vanessa = new TH1D("h_e_d_vanessa", "Delayed energy (Vanessa)", e_d_bins.size() - 1, e_d_bins.data());

        double e_dt_min = 0.0;
        double e_dt_max = 1.0;
        double e_dt_width = 0.025;
        int e_dt_nbin = std::round((e_dt_max - e_dt_min) / e_dt_width) + 1;
        std::vector<double> e_dt_bins = np::linspace(e_dt_min, e_dt_max, e_dt_nbin);
        TH1D* h_dt = new TH1D("h_dt", "Prompt-Delayed time difference", e_dt_bins.size() - 1, e_dt_bins.data());
        TH1D* h_dt_vanessa = new TH1D("h_dt_vanessa", "Prompt-Delayed time difference (Vanessa)", e_dt_bins.size() - 1, e_dt_bins.data());

        double e_dr_min = 0.0;
        double e_dr_max = 1.5;
        double e_dr_width = 0.05;
        int e_dr_nbin = std::round((e_dr_max - e_dr_min) / e_dr_width) + 1;
        std::vector<double> e_dr_bins = np::linspace(e_dr_min, e_dr_max, e_dr_nbin);
        TH1D* h_dr = new TH1D("h_dr", "Prompt-Delayed distance", e_dr_bins.size() - 1, e_dr_bins.data());

        double rho_min = 0.0;
        double rho_max = 17.7 * 17.7;
        int rho_nbin = 51;
        double z_min = -20.0;
        double z_max = 20.0;
        int z_nbin = 51;
        std::vector<double> rho_bins = np::linspace(rho_min, rho_max, rho_nbin);
        std::vector<double> z_bins = np::linspace(z_min, z_max, z_nbin);
        TH2D* h_rho_z_p = new TH2D("h_rho_z_p", "Prompt vertex distribution", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
        TH2D* h_rho_z_d = new TH2D("h_rho_z_d", "Delayed vertex distribution", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());

        for (const ibd& v : ibds) {
            h_e_p->Fill(v.prompt.e);
            h_e_d->Fill(v.delayed.e);
            h_dt->Fill(timestamp_to_double(v.delayed.ts - v.prompt.ts) * 1000.0);
            h_dr->Fill(mag(v.delayed.pos - v.prompt.pos) / 1000.0);
            h_rho_z_p->Fill((v.prompt.pos.x * v.prompt.pos.x + v.prompt.pos.y * v.prompt.pos.y) / 1.0e6, v.prompt.pos.z / 1000.0);
            h_rho_z_d->Fill((v.delayed.pos.x * v.delayed.pos.x + v.delayed.pos.y * v.delayed.pos.y) / 1.0e6, v.delayed.pos.z / 1000.0);
        }
        
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

        h_e_p_vanessa->SetLineStyle(3);
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

        h_e_d_vanessa->SetLineStyle(3);
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
    }

protected:

    std::shared_ptr<basic_navigator> m_nav;

    std::set<ibd> m_ibds;

};

using ibd_analysis = basic_analysis;

class cosmo_shape_analysis : public basic_analysis {

    // Essayer de voir si on peut, par le calcul, retrouver le nombre théorique de cosmo en prenant en compte les inéfficacité
    // Essayer de faire la même chose mais avec les autres méthodes de reconstruction, i.e. CdClassify et WpClassify

public:

    cosmo_shape_analysis(
        const std::string& filepath, const std::string& suffix, 
        const std::string& recname,
        const timestamp& sig_low, const timestamp& sig_high, 
        const timestamp& bkg_low, const timestamp& bkg_high, 
        double radius
    ) :
        basic_analysis{filepath, suffix},
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

        // double e_p_min = 0.7;
        // double e_p_max = 12.0;
        // double e_p_width = 0.20;
        // int e_p_nbin = std::round((e_p_max - e_p_min) / e_p_width) + 1;
        // std::vector<double> e_p_bins = np::linspace(e_p_min, e_p_max, e_p_nbin);
        std::vector<double> e_p_bins = create_custom_e_p_bins();
        TH1D* h_e_p_cosmo_bkg = new TH1D("h_e_p_cosmo_bkg", "Prompt energy (Cosmo bkg)", e_p_bins.size() - 1, e_p_bins.data());
        TH1D* h_e_p_cosmo_sig = new TH1D("h_e_p_cosmo_sig", "Prompt energy (Cosmo sig)", e_p_bins.size() - 1, e_p_bins.data());
        TH1D* h_e_p_cosmo_diff = new TH1D("h_e_p_cosmo_diff", "Prompt energy (Cosmo sig - Cosmo bkg)", e_p_bins.size() - 1, e_p_bins.data());

        double e_d_min = 2.0;
        double e_d_max = 2.5;
        double e_d_width = 0.02;
        int e_d_nbin = std::round((e_d_max - e_d_min) / e_d_width) + 1;
        std::vector<double> e_d_bins = np::linspace(e_d_min, e_d_max, e_d_nbin);
        TH1D* h_e_d_cosmo_bkg = new TH1D("h_e_d_cosmo_bkg", "Delayed energy (Cosmo bkg)", e_d_bins.size() - 1, e_d_bins.data());
        TH1D* h_e_d_cosmo_sig = new TH1D("h_e_d_cosmo_sig", "Delayed energy (Cosmo sig)", e_d_bins.size() - 1, e_d_bins.data());

        double e_dt_min = 0.0;
        double e_dt_max = 1.0;
        double e_dt_width = 0.025;
        int e_dt_nbin = std::round((e_dt_max - e_dt_min) / e_dt_width) + 1;
        std::vector<double> e_dt_bins = np::linspace(e_dt_min, e_dt_max, e_dt_nbin);
        TH1D* h_dt_cosmo_bkg = new TH1D("h_dt_cosmo_bkg", "Prompt-Delayed time difference (Cosmo bkg)", e_dt_bins.size() - 1, e_dt_bins.data());
        TH1D* h_dt_cosmo_sig = new TH1D("h_dt_cosmo_sig", "Prompt-Delayed time difference (Cosmo sig)", e_dt_bins.size() - 1, e_dt_bins.data());

        double e_dr_min = 0.0;
        double e_dr_max = 1.5;
        double e_dr_width = 0.05;
        int e_dr_nbin = std::round((e_dr_max - e_dr_min) / e_dr_width) + 1;
        std::vector<double> e_dr_bins = np::linspace(e_dr_min, e_dr_max, e_dr_nbin);
        TH1D* h_dr_cosmo_bkg = new TH1D("h_dr_cosmo_bkg", "Prompt-Delayed distance (Cosmo bkg)", e_dr_bins.size() - 1, e_dr_bins.data());
        TH1D* h_dr_cosmo_sig = new TH1D("h_dr_cosmo_sig", "Prompt-Delayed distance (Cosmo sig)", e_dr_bins.size() - 1, e_dr_bins.data());

        double rho_min = 0.0;
        double rho_max = 17.7 * 17.7;
        int rho_nbin = 51;
        double z_min = -20.0;
        double z_max = 20.0;
        int z_nbin = 51;
        std::vector<double> rho_bins = np::linspace(rho_min, rho_max, rho_nbin);
        std::vector<double> z_bins = np::linspace(z_min, z_max, z_nbin);
        TH2D* h_rho_z_p_cosmo_bkg = new TH2D("h_rho_z_p_cosmo_bkg", "Prompt vertex distribution (Cosmo bkg)", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
        TH2D* h_rho_z_d_cosmo_bkg = new TH2D("h_rho_z_d_cosmo_bkg", "Delayed vertex distribution (Cosmo bkg)", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
        TH2D* h_rho_z_p_cosmo_sig = new TH2D("h_rho_z_p_cosmo_sig", "Prompt vertex distribution (Cosmo sig)", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
        TH2D* h_rho_z_d_cosmo_sig = new TH2D("h_rho_z_d_cosmo_sig", "Delayed vertex distribution (Cosmo sig)", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());

        for (const cosmogenic& v : cosmos_bkg) {
            h_e_p_cosmo_bkg->Fill(v.prompt.e);
            h_e_d_cosmo_bkg->Fill(v.delayed.e);
            h_dt_cosmo_bkg->Fill(timestamp_to_double(v.delayed.ts - v.prompt.ts) * 1000.0);
            h_dr_cosmo_bkg->Fill(mag(v.delayed.pos - v.prompt.pos) / 1000.0);
            h_rho_z_p_cosmo_bkg->Fill((v.prompt.pos.x * v.prompt.pos.x + v.prompt.pos.y * v.prompt.pos.y) / 1.0e6, v.prompt.pos.z / 1000.0);
            h_rho_z_d_cosmo_bkg->Fill((v.delayed.pos.x * v.delayed.pos.x + v.delayed.pos.y * v.delayed.pos.y) / 1.0e6, v.delayed.pos.z / 1000.0);
        }

        for (const cosmogenic& v : cosmos_sig) {
            h_e_p_cosmo_sig->Fill(v.prompt.e);
            h_e_d_cosmo_sig->Fill(v.delayed.e);
            h_dt_cosmo_sig->Fill(timestamp_to_double(v.delayed.ts - v.prompt.ts) * 1000.0);
            h_dr_cosmo_sig->Fill(mag(v.delayed.pos - v.prompt.pos) / 1000.0);
            h_rho_z_p_cosmo_sig->Fill((v.prompt.pos.x * v.prompt.pos.x + v.prompt.pos.y * v.prompt.pos.y) / 1.0e6, v.prompt.pos.z / 1000.0);
            h_rho_z_d_cosmo_sig->Fill((v.delayed.pos.x * v.delayed.pos.x + v.delayed.pos.y * v.delayed.pos.y) / 1.0e6, v.delayed.pos.z / 1000.0);
        }

        h_e_p_cosmo_diff->Add(h_e_p_cosmo_sig, h_e_p_cosmo_bkg, 1.0, -1.0);

        // ============================================================================================
        // Cosmo bkg - Prompt energy
        // ============================================================================================

        TCanvas* c_e_p_cosmo_bkg = new TCanvas("c_e_p_cosmo_bkg", "Prompt energy (Cosmo bkg)", 1000, 1000);
        c_e_p_cosmo_bkg->cd();

        h_e_p_cosmo_bkg->SetLineWidth(3);
        h_e_p_cosmo_bkg->SetLineStyle(kSolid);
        h_e_p_cosmo_bkg->SetLineColorAlpha(kBlue, 1.0);

        h_e_p_cosmo_bkg->Draw();

        c_e_p_cosmo_bkg->Update();

        // ============================================================================================
        // Cosmo bkg - Delayed energy
        // ============================================================================================

        TCanvas* c_e_d_cosmo_bkg = new TCanvas("c_e_d_cosmo_bkg", "Delayed energy (Cosmo bkg)", 1000, 1000);
        c_e_d_cosmo_bkg->cd();

        h_e_d_cosmo_bkg->SetLineWidth(3);
        h_e_d_cosmo_bkg->SetLineStyle(kSolid);
        h_e_d_cosmo_bkg->SetLineColorAlpha(kBlue, 1.0);

        h_e_d_cosmo_bkg->Draw();

        c_e_d_cosmo_bkg->Update();

        // ============================================================================================
        // Cosmo bkg - Prompt-Delayed time difference
        // ============================================================================================

        TCanvas* c_dt_cosmo_bkg = new TCanvas("c_dt_cosmo_bkg", "Prompt-Delayed time difference (Cosmo bkg)", 1000, 1000);
        c_dt_cosmo_bkg->cd();
    
        h_dt_cosmo_bkg->SetLineWidth(3);
        h_dt_cosmo_bkg->SetLineStyle(kSolid);
        h_dt_cosmo_bkg->SetLineColorAlpha(kBlue, 1.0);

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

class cosmo_rate_analysis : public basic_analysis {

    // Regarder nombre cosmo en fonction du nombre de neutron ==> est-ce que ça suit une loi de poisson

public:

    cosmo_rate_analysis(const std::string& filepath, const std::string& suffix) :
        basic_analysis{filepath, suffix}
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
        TH1D* h_cosmo_rate_with_neutron = new TH1D("h_cosmo_rate_with_neutron", "Cosmo Rate With Neutron", 120, 0.0, 1.2);
        TH1D* h_cosmo_rate_with_no_neutron = new TH1D("h_cosmo_rate_with_no_neutron", "Cosmo Rate With No Neutron", 120, 0.0, 1.2);
        TH1D* h_cosmo_rate_with_at_least_1_neutron = new TH1D("h_cosmo_rate_with_at_least_1_neutron", "Cosmo Rate With At Least 1 Neutron", 120, 0.0, 1.2);
        TH1D* h_cosmo_rate_with_at_least_2_neutron = new TH1D("h_cosmo_rate_with_at_least_2_neutron", "Cosmo Rate With At Least 2 Neutron", 120, 0.0, 1.2);
        TH1D* h_cosmo_rate_with_at_least_3_neutron = new TH1D("h_cosmo_rate_with_at_least_3_neutron", "Cosmo Rate With At Least 3 Neutron", 120, 0.0, 1.2);
        for (const std::pair<ibd, std::vector<muon_data_association>>& val : m_ibds_to_mu) {
            const std::vector<muon_data_association>& muon_data = val.second;
            for (const muon_data_association& assoc : muon_data) {
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

        // TF1* f = new TF1(Form("f_%s", h->GetName()), "[0] + [1] * exp(-x / [2])", 0.02, 1.2);
        // f->SetParameter(0, constant_term);
        // f->SetParameter(1, exponential_term);
        // f->SetParameter(2, 180.0e-3);

        TF1* f = new TF1(Form("f_%s", h->GetName()), "[0] * [1] * [2] * exp(- [2] * x) + (1 - [0]) * [1] * [3] * exp(- [3] * x) + [4]", 0.02, 1.2);
        // f * N * lambda_Li9 * exp(- lambda_Li9 * t) + (1 - f) * N * lamnda_He8 * exp(- lambda_He8 * t) + c
        f->SetParameter(0, 0.85);
        f->SetParameter(1, exponential_term);
        f->SetParameter(2, 178.0e-3);
        f->SetParameter(3, 119.0e-3);
        f->SetParameter(4, constant_term);

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
        std::cout << "nIBD = " << f->GetParameter(4) * time_window / binning << " +/- " << f->GetParError(4) * time_window / binning << '\n';
        std::cout << "nLi = " << f->GetParameter(0) * f->GetParameter(1) << '\n';
        std::cout << "nHe = " << (1.0 - f->GetParameter(0)) * f->GetParameter(1) << '\n';
        std::cout << "Li decay = " << f->GetParameter(2) << '\n';
        std::cout << "He decay = " << f->GetParameter(3) << '\n';
        std::cout << "Fit Results for " << h->GetName() << ":\n";


        // double time_window = 1.2;
        // double binning = time_window / 120.0;
        // std::cout << "Fit Results for " << h->GetName() << ":\n";
        // std::cout << "nIBD = " << f->GetParameter(0) * time_window / binning << " +/- " << f->GetParError(0) * time_window / binning << '\n';
        // std::cout << "nLiHe = " << f->GetParameter(1) * f->GetParameter(2) * (1 - std::exp(-time_window / f->GetParameter(2))) / binning << '\n'; 


    }

};

class analysis_registry {

public:

    typedef std::vector<std::shared_ptr<analysis_base>> vector_type;
    typedef std::map<std::shared_ptr<navigator_base>, vector_type> map_type;

    bool book(const std::shared_ptr<analysis_base>& analysis) {
        if (!analysis) {
            std::cerr << "Cannot register analysis\n";
            return false;
        }
        m_registry[analysis->navigator()].push_back(analysis);
        return true;
    }

private:

    map_type m_registry;

    friend class analysis_manager;

};

class analysis_manager {

public:

    analysis_manager(analysis_registry& reg) :
        m_reg{reg}
    {}

    bool run() {
        if (m_reg.m_registry.empty()) {
            std::cout << "Analysis registry is empty. Exiting run\n";
            return false;
        }
        
        std::cout << "\n--- Starting Analysis Loop Over " << m_reg.m_registry.size() << " Data Groups ---\n";
        for (auto const& [nav, analyses] : m_reg.m_registry) {
            if (!nav->is_valid()) {
                std::cerr << "\nWARNING: Navigator for this group is invalid. Skipping group\n";
                continue;
            }

            Long64_t entries = nav->size();
            std::cout << "\n[Group Start] Processing " << analyses.size() << " analyses over " << entries << " entries\n";
            
            for (Long64_t k = 0; k < entries; ++k) {
                nav->entry(k); 

                for (const auto& analysis : analyses) {
                    if (analysis->selection()) {
                        if (!analysis->process()) return false;
                    }
                }

                if ((k + 1) % 1000 == 0) {
                     std::cout << "  [Group Status] Processed entry " << k + 1 << " / " << entries << '\n';
                }
            }

            for (const auto& analysis : analyses) {
                analysis->result();
            }
            std::cout << "[Group End] Finished processing group\n";
        }
        std::cout << "\n--- All Analysis Groups Finished ---\n";
        return true;
    }

private:

    analysis_registry& m_reg;

};

#endif // ANALYSIS_HPP_