#ifndef ANALYSIS_COSMO_SHAPE_ANALYSIS_HPP_
#define ANALYSIS_COSMO_SHAPE_ANALYSIS_HPP_

#include <TLegend.h>

#include "analysis/basic_analysis.hpp"
#include "utils/muon_lookup.hpp"
#include "utils/plot.hpp"

class cosmo_shape_analysis : public basic_analysis {

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
        double e_p = m_nav->prompt.e; // / m_gtc.interpolate(m_nav->prompt.ts);
        double e_d = m_nav->delayed.e; // / m_gtc.interpolate(m_nav->delayed.ts);

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
            double e_mult = m_nav->e_mult[k]; // / m_gtc.interpolate(ts_mult);
            if (e_mult < 2.0 || 12.0 < e_mult) continue;
            if (ts_mult < m_nav->prompt.ts - timestamp{0, 1000000} || m_nav->delayed.ts + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }
        if (nb_multu_veto) return false;

        if ( std::pow((m_nav->meta_prompt.stdhit - 0.55) / 0.45, 2.0) + std::pow((m_nav->meta_prompt.stdt - 170.0) / 80.0, 2.0) > 1.0 ) return false;

        m_dlat_mu2p.clear();
        m_dlat_mu2d.clear();
        m_dt_mu2p.clear();
        m_dt_mu2d.clear();
        m_is_sig.clear();

        multiplicity_muon_lookup nb_muons_in_cd_event;
        multiplicity_muon_lookup nb_muons_in_wp_event;
        nb_muons_in_cd_event.fill(m_nav, "CdClassify");
        nb_muons_in_wp_event.fill(m_nav, "WpBasic");

        // stopping_muon_lookup has_stopping_in_cd_event;
        // stopping_muon_lookup has_stopping_in_wp_event;
        // has_stopping_in_cd_event.fill(m_nav, "CdClassify");
        // has_stopping_in_wp_event.fill(m_nav, "WpBasic");

        for (std::size_t k = 0ul; k < m_nav->method_mu.size(); ++k) {
            if (m_nav->method_mu[k] != m_recname) continue;
            timestamp ts_mu{m_nav->sec_mu[k], m_nav->nsec_mu[k]};
            if (nb_muons_in_cd_event[ts_mu] > 1ul || nb_muons_in_wp_event[ts_mu] > 1ul) continue;
            // if (has_stopping_in_cd_event[ts_mu]) continue;
            // if (has_stopping_in_wp_event[ts_mu]) continue;
            
            bool found_neutron = false;
            for (std::size_t l = 0ul; l < m_nav->e_n.size() && !found_neutron; ++l) {
                timestamp ts_n{m_nav->sec_n[l], m_nav->nsec_n[l]};
                if (ts_n < ts_mu + timestamp{0, 20000} || ts_mu + timestamp{0, 2000000} < ts_n) continue;
                found_neutron = true;
            }
            if (!found_neutron) continue;

            bool is_in_bkg = (
                ts_mu + m_ts_bkg_low <= m_nav->prompt.ts && m_nav->prompt.ts <= ts_mu + m_ts_bkg_high &&
                ts_mu + m_ts_bkg_low <= m_nav->delayed.ts && m_nav->delayed.ts <= ts_mu + m_ts_bkg_high
            );
            bool is_in_sig = (
                ts_mu + m_ts_sig_low <= m_nav->prompt.ts && m_nav->prompt.ts <= ts_mu + m_ts_sig_high &&
                ts_mu + m_ts_sig_low <= m_nav->delayed.ts && m_nav->delayed.ts <= ts_mu + m_ts_sig_high
            );
            if (!is_in_bkg && !is_in_sig) continue;
            vec3 pos_mu{m_nav->posx_mu[k], m_nav->posy_mu[k], m_nav->posz_mu[k]};
            vec3 dir_mu = unit(vec3{m_nav->dirx_mu[k], m_nav->diry_mu[k], m_nav->dirz_mu[k]});
            double d_mu2p = mag(cross(dir_mu, m_nav->prompt.pos - pos_mu));
            double d_mu2d = mag(cross(dir_mu, m_nav->delayed.pos - pos_mu));
            if (m_radius < d_mu2p && m_radius < d_mu2d) continue;
            m_dlat_mu2p.push_back(d_mu2p);
            m_dlat_mu2d.push_back(d_mu2d);
            m_dt_mu2p.push_back(timestamp_to_double(m_nav->prompt.ts - ts_mu));
            m_dt_mu2d.push_back(timestamp_to_double(m_nav->delayed.ts - ts_mu));
            m_is_sig.push_back(is_in_sig);
        }

        return !m_is_sig.empty();
    }

    bool process() override {
        for (std::size_t k = 0ul; k < m_is_sig.size(); ++k) {
            if (m_is_sig[k]) {
                m_cosmos_sig.insert({m_nav->run_id, m_nav->prompt, m_nav->delayed, m_dlat_mu2p[k], m_dlat_mu2d[k], m_dt_mu2p[k], m_dt_mu2d[k]});
            }
            else {
                m_cosmos_bkg.insert({m_nav->run_id, m_nav->prompt, m_nav->delayed, m_dlat_mu2p[k], m_dlat_mu2d[k], m_dt_mu2p[k], m_dt_mu2d[k]});
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
            e_p = it->prompt.e; // / m_gtc.interpolate(it->prompt.ts);
            e_d = it->delayed.e; // / m_gtc.interpolate(it->delayed.ts);
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
            e_p = it->prompt.e; // / m_gtc.interpolate(it->prompt.ts);
            e_d = it->delayed.e; // / m_gtc.interpolate(it->delayed.ts);
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

    void result() override {}

protected:

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