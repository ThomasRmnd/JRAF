#include <set>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include "analysis/analysis.hpp"
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

struct VanessaIBD {

    int run_id;

    timestamp ts_p;
    double e_p;
    double totq_p;

    timestamp ts_d;
    double e_d;
    double totq_d;

};

inline bool operator<(const VanessaIBD& lhs, const VanessaIBD& rhs) {
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