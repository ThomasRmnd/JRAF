#ifndef ANALYSIS_IBD_ANALYSIS_HPP_
#define ANALYSIS_IBD_ANALYSIS_HPP_

#include <set>

#include <TFile.h>
#include <TTree.h>

#include "analysis/basic_analysis.hpp"
#include "utils/muon_lookup.hpp"

struct ibd_wmu {

    ibd i;
    timestamp dt_last_mu;
    double dlat_mu2p; // to closest muon
    timestamp dt_mu2p; // to closest muon

};

inline bool operator<(const ibd_wmu& lhs, const ibd_wmu& rhs) {
    return lhs.i < rhs.i;
}

class ibd_analysis : public basic_analysis {

public:

    ibd_analysis(const std::string& name, const std::string& filepath, const std::string& suffix) :
        basic_analysis{name, filepath, suffix}
    {}

    virtual ~ibd_analysis() override = default;

    virtual bool selection() override = 0;

    bool process() override {
        m_ibds.insert({
            {m_nav->run_id, m_nav->prompt, m_nav->delayed}, 
            m_dt_last_mu, m_dlat_mu2p, m_dt_mu2p
        });
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
        double dt_last_mu;
        double dlat_mu2p;
        double dt_mu2p;
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
        t->Branch("dt_last_mu", &dt_last_mu);
        t->Branch("dlat_mu2p", &dlat_mu2p);
        t->Branch("dt_mu2p", &dt_mu2p);

        for (std::set<ibd_wmu>::const_iterator it = m_ibds.begin(); it != m_ibds.end(); ++it) {
            run_id = it->i.run_id;
            pos_p = it->i.prompt.pos;
            pos_d = it->i.delayed.pos;
            ts_p = it->i.prompt.ts;
            ts_d = it->i.delayed.ts;
            e_p = it->i.prompt.e / m_gtc.interpolate(it->i.prompt.ts);
            e_d = it->i.delayed.e / m_gtc.interpolate(it->i.delayed.ts);
            dt_last_mu = timestamp_to_double(it->dt_last_mu);
            dlat_mu2p = it->dlat_mu2p;
            dt_mu2p = timestamp_to_double(it->dt_mu2p);
            t->Fill();
        }

        f->cd();
        t->Write();
        f->Close();
        return true;
    }

    void result() override {}

protected:

    std::set<ibd_wmu> m_ibds;
    timestamp m_dt_last_mu;
    double m_dlat_mu2p;
    timestamp m_dt_mu2p;

    void calculate_dt_to_last_muon() {
        m_dt_last_mu = timestamp{0, 0};
        bool is_set_dt_last_mu = false;
        for (std::size_t k = 0ul; k < m_nav->method_mu.size(); ++k) {
            timestamp ts_mu{m_nav->sec_mu[k], m_nav->nsec_mu[k]};
            if (m_nav->prompt.ts < ts_mu) continue;
            bool found_neutron = false;
            for (std::size_t l = 0ul; l < m_nav->e_n.size() && !found_neutron; ++l) {
                timestamp ts_n{m_nav->sec_n[l], m_nav->nsec_n[l]};
                if (ts_n < ts_mu + timestamp{0, 20000} || ts_mu + timestamp{0, 2000000} < ts_n) continue;
                found_neutron = true;
            }
            if (!found_neutron) continue;
            if (is_set_dt_last_mu && m_nav->prompt.ts - ts_mu > m_dt_last_mu) continue;
            m_dt_last_mu = m_nav->prompt.ts - ts_mu;
            is_set_dt_last_mu = true;
        }
        if (!is_set_dt_last_mu) {
            m_dt_last_mu = timestamp{-1, 0};
        }
    }

    void calculate_dlat_dt_muon_to_prompt() {
        multiplicity_muon_lookup nb_muons_in_cd_event;
        multiplicity_muon_lookup nb_muons_in_wp_event;
        nb_muons_in_cd_event.fill(m_nav, "CdClassify");
        nb_muons_in_wp_event.fill(m_nav, "WpBasic");

        // stopping_muon_lookup has_stopping_in_cd_event;
        // stopping_muon_lookup has_stopping_in_wp_event;
        // has_stopping_in_cd_event.fill(m_nav, "CdClassify");
        // has_stopping_in_wp_event.fill(m_nav, "WpBasic");

        m_dlat_mu2p = std::numeric_limits<double>::infinity();
        m_dt_mu2p = timestamp{0, 0};
        bool is_set_dlat_mu2p = false;
        for (std::size_t k = 0ul; k < m_nav->method_mu.size(); ++k) {
            if (m_nav->method_mu[k] != "CdWpTtChi2") continue;
            timestamp ts_mu{m_nav->sec_mu[k], m_nav->nsec_mu[k]};
            if (nb_muons_in_cd_event[ts_mu] > 1ul || nb_muons_in_wp_event[ts_mu] > 1ul) continue;
            // if (has_stopping_in_cd_event[ts_mu]) continue;
            // if (has_stopping_in_wp_event[ts_mu]) continue;

            bool is_in_sig = (
                ts_mu + timestamp{0, 5000000} <= m_nav->prompt.ts && m_nav->prompt.ts <= ts_mu + timestamp{0, 1200000000} &&
                ts_mu + timestamp{0, 5000000} <= m_nav->delayed.ts && m_nav->delayed.ts <= ts_mu + timestamp{0, 1200000000}
            );

            if (!is_in_sig) continue;
#ifdef __USE_MUON_IPOS_FPOS__
            vec3 pos_mu{m_nav->iposx_mu[k], m_nav->iposy_mu[k], m_nav->iposz_mu[k]};
            vec3 fpos_mu{m_nav->fposx_mu[k], m_nav->fposy_mu[k], m_nav->fposz_mu[k]};
            vec3 dir_mu = unit(fpos_mu - pos_mu);
#else
            vec3 pos_mu{m_nav->posx_mu[k], m_nav->posy_mu[k], m_nav->posz_mu[k]};
            vec3 dir_mu = unit(vec3{m_nav->dirx_mu[k], m_nav->diry_mu[k], m_nav->dirz_mu[k]});
#endif
            if (
                std::isnan(pos_mu.x) || std::isnan(pos_mu.y) || std::isnan(pos_mu.z) ||
                std::isnan(dir_mu.x) || std::isnan(dir_mu.y) || std::isnan(dir_mu.z)
            ) continue;

            double d_mu2p = mag(cross(dir_mu, m_nav->prompt.pos - pos_mu));

            if (m_dlat_mu2p < d_mu2p) continue;
            m_dlat_mu2p = d_mu2p;
            m_dt_mu2p = m_nav->prompt.ts - ts_mu;
            is_set_dlat_mu2p = true;
        }

        if (!is_set_dlat_mu2p) {
            m_dlat_mu2p = -1.0;
            m_dt_mu2p = timestamp{-1, 0};
        }
    }

};

#endif // ANALYSIS_IBD_ANALYSIS_HPP_