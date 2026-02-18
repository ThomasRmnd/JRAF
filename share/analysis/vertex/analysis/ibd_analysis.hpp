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

struct ibd_wmu {

    ibd i;
    double dt_last_mu;

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

    virtual bool selection() override {
        double e_p = m_nav->prompt.e / m_gtc.interpolate(m_nav->prompt.ts);
        double e_d = m_nav->delayed.e / m_gtc.interpolate(m_nav->delayed.ts);

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
            double e_mult = m_nav->e_mult[k] / m_gtc.interpolate(ts_mult);
            if (e_mult < 2.0 || 12.0 < e_mult) continue;
            if (ts_mult < m_nav->prompt.ts - timestamp{0, 1000000} || m_nav->delayed.ts + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }
        if (nb_multu_veto) return false;

        std::size_t nb_neutron_veto = 0ul;
        for (std::size_t k = 0ul; k < m_nav->e_n.size(); ++k) {
            timestamp ts_n{m_nav->sec_n[k], m_nav->nsec_n[k]};
            vec3 pos_n{m_nav->posx_n[k], m_nav->posy_n[k], m_nav->posz_n[k]};
            double e_n = m_nav->e_n[k] / m_gtc.interpolate(ts_n);
            if (e_n < 1.5 || 20.0 < e_n) continue;
            if (m_nav->stdt_n[k] > 275.0) continue;
            if (
                (mag(m_nav->prompt.pos - pos_n) < 4000.0 && ts_n + timestamp{0, 20000} < m_nav->prompt.ts && m_nav->prompt.ts < ts_n + timestamp{0, 1200000000}) ||
                (mag(m_nav->delayed.pos - pos_n) < 4000.0 && ts_n + timestamp{0, 20000} < m_nav->delayed.ts && m_nav->delayed.ts < ts_n + timestamp{0, 1200000000})
            ) {
                ++nb_neutron_veto;
            }
        }
        if (nb_neutron_veto) return false;

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
            if (!is_set_dt_last_mu && m_nav->prompt.ts - ts_mu > m_dt_last_mu) continue;
            m_dt_last_mu = m_nav->prompt.ts - ts_mu;
            is_set_dt_last_mu = true;
        }

        if ( std::pow((m_nav->meta_prompt.stdhit - 0.55) / 0.45, 2.0) + std::pow((m_nav->meta_prompt.stdt - 170.0) / 80.0, 2.0) > 1.0 ) return false;

        return true;
    }

    bool process() override {
        m_ibds.insert({{m_nav->run_id, m_nav->prompt, m_nav->delayed}, timestamp_to_double(m_dt_last_mu)});
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

        for (std::set<ibd_wmu>::const_iterator it = m_ibds.begin(); it != m_ibds.end(); ++it) {
            run_id = it->i.run_id;
            pos_p = it->i.prompt.pos;
            pos_d = it->i.delayed.pos;
            ts_p = it->i.prompt.ts;
            ts_d = it->i.delayed.ts;
            e_p = it->i.prompt.e / m_gtc.interpolate(it->i.prompt.ts);
            e_d = it->i.delayed.e / m_gtc.interpolate(it->i.delayed.ts);
            dt_last_mu = it->dt_last_mu;
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

};

#endif // ANALYSIS_IBD_ANALYSIS_HPP_