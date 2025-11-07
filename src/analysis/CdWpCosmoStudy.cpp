#include "analysis/CdWpCosmoStudy.hpp"

#include <algorithm>

#include "SniperKernel/SniperLog.h"

#include "event/IBD.hpp"
#include "selection/Energy.hpp"
#include "selection/Muon.hpp"
#include "selection/Vertex.hpp"
#include "selection/Volume.hpp"

CdWpCosmoStudy::CdWpCosmoStudy(const std::string& name, const std::string& method, double cyl_radius, const TimeStamp& lwr_window, const TimeStamp& upr_window) : 
    Analysis{name, method}, 
    m_cyl_radius{cyl_radius},
    m_lwr_window{lwr_window}, 
    m_upr_window{upr_window} 
{}

bool CdWpCosmoStudy::initialize() {
    if (!Analysis::initialize()) return false;
    m_tree->Branch("dlat_p", &dlat_p);
    m_tree->Branch("dlat_d", &dlat_d);
    m_tree->Branch("dt_mu2p_sec", &dt_mu2p_sec);
    m_tree->Branch("dt_mu2d_sec", &dt_mu2d_sec);
    m_tree->Branch("dt_mu2p_nsec", &dt_mu2p_nsec);
    m_tree->Branch("dt_mu2d_nsec", &dt_mu2d_nsec);
    return true; 
}

void CdWpCosmoStudy::process(const EventContext::View& events) {
    std::vector<TimeRangeMuonVetoSelection> mu_wp_bundle_cut;
    std::vector<CylindricalMuonVetoSelection> mu_cosmo_cut;
    for (const track& trk : events.tracks()) {
        std::vector<track> cd_tracks;
        std::vector<track> wp_tracks;
        std::vector<track> tt_tracks;
        bool has_cdclassify = false;
        bool has_cdwpttchi2 = false;
        for (const track& trk2 : events.tracks()) {
            mu_wp_bundle_cut.emplace_back(trk2, TimeStamp{0, 0}, TimeStamp{0, 5000000});
            if (trk2.method == "CdClassify") has_cdclassify = true;
            if (trk2.method == "CdWpTtChi2") has_cdwpttchi2 = true;
            if (trk2.method != "WpBasic" /* jt->det != track::loc::wp */) continue;
            wp_tracks.push_back(trk2);
        }
        bool has_cdclassify_but_no_cdwpttchi2 = has_cdclassify && !has_cdwpttchi2;
        if (wp_tracks.size() > 1ul) {
            mu_wp_bundle_cut.emplace_back(trk, TimeStamp{0, 0}, TimeStamp{0, 500000000});
            // std::cout << "[DEBUG] Bundle veto from starting at " << it->begin()->ts << '\n';
            continue;
        }
        for (const track& trk2 : events.tracks()) {
            if (trk2.method != "CdWpTtChi2" /* jt->det != track::loc::cd */) continue;
            if (trk2.quality != -1.0f && !has_cdclassify_but_no_cdwpttchi2) {
                cd_tracks.push_back(trk2);
                mu_cosmo_cut.emplace_back(trk2, m_cyl_radius, m_lwr_window, m_upr_window);
            }
            else {
                mu_wp_bundle_cut.emplace_back(trk2, TimeStamp{0, 0}, TimeStamp{0, 500000000});
            }
            // else if (jt->det == track::loc::tt) {
            //     tt_tracks.push_back(*jt);
            // }
        }
    }

    FiducialVolumeSelection fiducial_vol_cut{16500.0};
    ChimneySelection chimney_cut{15500.0, 3000.0};
    // ChargeRangeSelection prompt_charge_cut{1500.0, 20000.0};
    // ChargeRangeSelection delayed_charge_cut{4000.0, 6000.0};
    EnergyRangeSelection prompt_energy_cut{0.7, 12.0};
    EnergyRangeSelection delayed_energy_cut{2.0, 2.5};
    EnergyRangeSelection multiplicity_energy_cut{2.0, 12.0};

    std::vector<std::vector<track>> tracks_for_cosmo;
    std::vector<ibd> cosmos;

    for (const vertex& prompt : events.current()) {
        LogInfo << prompt << '\n';
        if (!fiducial_vol_cut.isIn(prompt)) {
            LogInfo << "Prompt not in fiducial volume\n";
            continue;
        }
        if (chimney_cut.isIn(prompt)) {
            LogInfo << "Prompt is a chimney\n";
            continue;
        }

        if (!prompt_energy_cut.isIn(prompt)) {
            LogInfo << "Prompt not in energy range\n";
            continue;
        }

        bool is_vetoed = false;
        for (const TimeRangeMuonVetoSelection& cut : mu_wp_bundle_cut) {
            if (!cut.isIn(prompt)) continue;
            is_vetoed = true;
            break;
        }
        if (is_vetoed) {
            LogInfo << "Prompt is muon vetoed\n";
            continue;
        }

        is_vetoed = false;
        std::vector<CylindricalMuonVetoSelection> trk_cut_cand;
        for (const CylindricalMuonVetoSelection& cut : mu_cosmo_cut) {
            if (!cut.isIn(prompt)) continue;
            is_vetoed = true;
            trk_cut_cand.push_back(cut);
        }
        if (!is_vetoed) {
            LogInfo << "Prompt is not in cylindrical muon cosmogenic cut\n";
            continue;
        }

        TimeRangeSelection multi_prompt_time{prompt.ts, TimeStamp{0, -1000000}, TimeStamp{0, 0}};
        bool prompt_has_multi = false;
        for (const vertex& cand : events.before()) {
            if (!multi_prompt_time.isIn(cand)) continue;
            // if (!fiducial_vol_cut.isIn(cand)) continue;
            // if (chimney_cut.isIn(cand)) continue;
            if (!multiplicity_energy_cut.isIn(cand)) continue;
            is_vetoed = false;
            for (const TimeRangeMuonVetoSelection& cut : mu_wp_bundle_cut) {
                if (!cut.isIn(cand)) continue;
                is_vetoed = true;
                break;
            }
            if (is_vetoed) continue;
            is_vetoed = false;
            for (const CylindricalMuonVetoSelection& cut : mu_cosmo_cut) {
                if (!cut.isIn(cand)) continue;
                is_vetoed = true;
                break;
            }
            if (is_vetoed) continue;

            prompt_has_multi = true;
            break;
        }
        if (prompt_has_multi) {
            LogInfo << "Prompt has multiplicity\n";
            continue;
        }

        VertexCorrelationSelection correlation_cut{prompt, 1500.0, TimeStamp{0, 5000}, TimeStamp{0, 1000000}};

        for (const vertex& delayed : events.after()) {
            LogInfo << delayed << '\n';
            if (!fiducial_vol_cut.isIn(delayed)) {
                LogInfo << "Delayed not in fiducial volume\n";
                continue;
            }
            // if (chimney_cut.isIn(delayed)) {
            //     LogInfo << "Delayed is a chimney\n";
            //     continue;
            // }

            if (!delayed_energy_cut.isIn(delayed)) {
                LogInfo << "Delayed not in energy range\n";
                continue;
            }

            if (!correlation_cut.isIn(delayed)) {
                LogInfo << "Delayed not correlated\n";
                continue;
            }

            is_vetoed = false;
            for (const TimeRangeMuonVetoSelection& cut : mu_wp_bundle_cut) {
                if (!cut.isIn(delayed)) continue;
                is_vetoed = true;
                break;
            }
            if (is_vetoed) {
                LogInfo << "Delayed is muon vetoed\n";
                continue;
            }

            is_vetoed = false;
            // for (const CylindricalMuonVetoSelection& cut : mu_cosmo_cut) {
            std::vector<track> trk_cand;
            for (const CylindricalMuonVetoSelection& cut : trk_cut_cand) {
                if (!cut.isIn(delayed)) continue;
                is_vetoed = true;
                trk_cand.push_back(cut.c_trk);
            }
            if (!is_vetoed) {
                LogInfo << "Delayed is not in cylindrical muon cosmogenic cut\n";
                continue;
            }

            TimeRangeSelection multi_delayed_time{delayed.ts, TimeStamp{0, 0}, TimeStamp{0, 1000000}};
            bool delayed_has_multi = false;
            for (const vertex& cand : events.after()) {
                if (cand.ts == delayed.ts) continue; // same event
                // if (!fiducial_vol_cut.isIn(cand)) continue;
                // if (chimney_cut.isIn(cand)) continue;
                if (!multiplicity_energy_cut.isIn(cand)) continue;
                is_vetoed = false;
                for (const TimeRangeMuonVetoSelection& cut : mu_wp_bundle_cut) {
                    if (!cut.isIn(cand)) continue;
                    is_vetoed = true;
                    break;
                }
                if (is_vetoed) continue;
                is_vetoed = false;
                for (const CylindricalMuonVetoSelection& cut : mu_cosmo_cut) {
                    if (!cut.isIn(cand)) continue;
                    is_vetoed = true;
                    break;
                }
                if (is_vetoed) continue;

                if (cand.ts < delayed.ts) {
                    delayed_has_multi = true;
                    LogInfo << "Delayed is in the in-between multiplicity cut by " << cand.ts << '\n';
                    break; // in-between p-d multiplicity
                }
                if (!multi_delayed_time.isIn(cand)) continue;
                delayed_has_multi = true;
                LogInfo << "Delayed is in the after multiplicity cut by " << cand.ts << '\n';
                break; // after p-d multiplicity
            }
            if (delayed_has_multi) {
                LogInfo << "Delayed has multiplicity\n";
                continue;
            }

            cosmos.emplace_back(prompt, delayed);
            tracks_for_cosmo.push_back(trk_cand);
            LogInfo << "IBD event detected!\n";
        }
    }

    for (std::size_t k = 0ul; k < cosmos.size(); ++k) {
        const ibd& ibd_ = cosmos[k];
        const std::vector<track>& trk_cand = tracks_for_cosmo[k];
        for (const track& trk_ : trk_cand) {
            posx_p = ibd_.prompt.pos.x;
            posy_p = ibd_.prompt.pos.y;
            posz_p = ibd_.prompt.pos.z;
            e_p = ibd_.prompt.energy;
            totq_p = ibd_.prompt.totq;
            sec_p = ibd_.prompt.ts.GetSec();
            nsec_p = ibd_.prompt.ts.GetNanoSec();
            posx_d = ibd_.delayed.pos.x;
            posy_d = ibd_.delayed.pos.y;
            posz_d = ibd_.delayed.pos.z;
            e_d = ibd_.delayed.energy;
            totq_d = ibd_.delayed.totq;
            sec_d = ibd_.delayed.ts.GetSec();
            nsec_d = ibd_.delayed.ts.GetNanoSec();

            vec3 trk_dir = unit(trk_.fpos - trk_.ipos);
            dlat_p = mag(cross(trk_dir, ibd_.prompt.pos - trk_.ipos));
            dlat_d = mag(cross(trk_dir, ibd_.delayed.pos - trk_.ipos));
            dt_mu2p_sec = (ibd_.prompt.ts - trk_.ts).GetSec();
            dt_mu2d_sec = (ibd_.delayed.ts - trk_.ts).GetSec();
            dt_mu2p_nsec = (ibd_.prompt.ts - trk_.ts).GetNanoSec();
            dt_mu2d_nsec = (ibd_.delayed.ts - trk_.ts).GetNanoSec();
            m_tree->Fill();
        }
    }
}