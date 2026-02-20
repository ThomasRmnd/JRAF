#ifndef ANALYSIS_COSMO_MUON_SHAPE_ANALYSIS_HPP_
#define ANALYSIS_COSMO_MUON_SHAPE_ANALYSIS_HPP_

#include <cmath>
#include <vector>
#include <set>

#include <TFile.h>
#include <TTree.h>

#include "analysis/basic_analysis.hpp"
#include "utils/muon_lookup.hpp"

class cosmo_shape_muon_analysis : public basic_analysis {

public:

    cosmo_shape_muon_analysis(
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

    virtual ~cosmo_shape_muon_analysis() override = default;

    bool selection() override = 0;

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
            e_p = it->prompt.e / m_gtc.interpolate(it->prompt.ts);
            e_d = it->delayed.e / m_gtc.interpolate(it->delayed.ts);
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
            e_p = it->prompt.e / m_gtc.interpolate(it->prompt.ts);
            e_d = it->delayed.e / m_gtc.interpolate(it->delayed.ts);
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

#endif // ANALYSIS_COSMO_MUON_SHAPE_ANALYSIS_HPP_