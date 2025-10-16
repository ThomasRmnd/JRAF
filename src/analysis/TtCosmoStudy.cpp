#include "analysis/TtCosmoStudy.hpp"

#include <algorithm>

#include "SniperKernel/SniperLog.h"

#include "event/Event.hpp"
#include "event/IBD.hpp"
#include "selection/Energy.hpp"
#include "selection/Muon.hpp"
#include "selection/Volume.hpp"

TtCosmoStudy::TtCosmoStudy(const std::string& name, const TimeStamp& lwr_window, const TimeStamp& upr_window) : 
    Analysis{name}, 
    m_lwr_window{lwr_window}, 
    m_upr_window{upr_window} 
{}

bool TtCosmoStudy::initialize() {
    if (!Analysis::initialize()) return false;
    m_tree->Branch("totq_p", &totq_p);
    m_tree->Branch("totq_d", &totq_d);

    m_tree->Branch("dlat_p", &dlat_p);
    m_tree->Branch("dlat_d", &dlat_d);
    m_tree->Branch("dt_mu2p_sec", &dt_mu2p_sec);
    m_tree->Branch("dt_mu2d_sec", &dt_mu2d_sec);
    m_tree->Branch("dt_mu2p_nsec", &dt_mu2p_nsec);
    m_tree->Branch("dt_mu2d_nsec", &dt_mu2d_nsec);
    return true; 
}

void TtCosmoStudy::process(JM::NavBuffer* buf) {
    std::vector<std::vector<track>> tracks;
    std::vector<vertex> cur_vertices;
    std::vector<vertex> bef_vertices;
    std::vector<vertex> aft_vertices;
    tracks.reserve(buf->size());
    bef_vertices.reserve(buf->size() / 2);
    aft_vertices.reserve(buf->size() / 2);
    for (JM::NavBuffer::Iterator it = buf->begin(); it != buf->end(); ++it) {
        JM::EvtNavigator* nav = it->get();
        Event evt;
        evt.load(nav);
        tracks.push_back(evt.tracks);
        if (it < buf->current()) {
            bef_vertices.insert(bef_vertices.end(), evt.vertices.begin(), evt.vertices.end());
        }
        else if (buf->current() < it) {
            aft_vertices.insert(aft_vertices.end(), evt.vertices.begin(), evt.vertices.end());
        }
        else {
            LogInfo << evt << '\n';
            cur_vertices.insert(cur_vertices.end(), evt.vertices.begin(), evt.vertices.end());
        }
    }

    std::vector<WaterPoolMuonVetoSelection> mu_wp_bundle_cut;
    std::vector<BasicMuonVetoSelection> mu_cosmo_cut;
    for (std::vector<std::vector<track>>::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
        if (it->empty()) continue;
        std::vector<track> cd_tracks;
        std::vector<track> wp_tracks;
        std::vector<track> tt_tracks;
        for (std::vector<track>::const_iterator jt = it->begin(); jt != it->end(); ++jt) {
            mu_wp_bundle_cut.emplace_back(*jt, TimeStamp{0, 5000000});
            if (jt->det == track::loc::cd) {
                if (jt->quality != -1.0f) cd_tracks.push_back(*jt);
            }
            else if (jt->det == track::loc::wp) {
                wp_tracks.push_back(*jt);
            }
            else if (jt->det == track::loc::tt) {
                tt_tracks.push_back(*jt);
            }
        }
        if (cd_tracks.size() != 1ul || wp_tracks.size() != 1ul || tt_tracks.size() != 1ul) continue;
        track cd_trk = cd_tracks.front();
        track tt_trk = tt_tracks.front();
        vec3 cd_dir = unit(cd_trk.fpos - cd_trk.ipos);
        vec3 tt_dir = unit(tt_trk.fpos - tt_trk.ipos);
        vec3 cd_mid = cd_trk.ipos - dot(cd_trk.ipos, cd_dir) * cd_dir;
        vec3 tt_mid = tt_trk.ipos - dot(tt_trk.ipos, tt_dir) * tt_dir;
        double angle_cdtt = angle(cd_dir, tt_dir) * 180 / M_PI;
        double distance_cdtt = mag(cd_mid - tt_mid);
        double clippingness = mag(cross(tt_dir, -tt_trk.ipos));

        if (
            (clippingness > 18.0 || angle_cdtt > 3.0 || distance_cdtt > 1.0) &&
            (clippingness < 16.0 || clippingness > 18.0 || angle_cdtt > 5.0 || distance_cdtt > 1.5)
        ) {
            continue;
        }
        mu_cosmo_cut.emplace_back(cd_trk, 3000.0, m_lwr_window, m_upr_window);
    }

    FiducialVolumeSelection fiducial_vol_cut{17200.0};
    HeightVolumeSelection lower_height_vol_cut{-20050.0, -11000.0};
    HeightVolumeSelection upper_height_vol_cut{ 11000.0,  20050.0};
    XYRadiusVolumeSelection xyradius_vol_cut{0.0, 3000.0};
    double prompt_lower_thold = 1500.0;
    double prompt_upper_thold = 20000.0;
    double delayed_lower_thold = 3700.0;
    double delayed_upper_thold = 6000.0;

    std::vector<ibd> cosmos;
    std::vector<track> tracks_for_cosmo;

    for (const vertex& prompt : cur_vertices) {
        LogInfo << prompt << '\n';
        if (!fiducial_vol_cut.isIn(prompt)) {
            LogInfo << "Prompt not in fiducial volume\n";
            continue;
        }
        
        if (
            (upper_height_vol_cut.isIn(prompt) && xyradius_vol_cut.isIn(prompt)) ||
            (lower_height_vol_cut.isIn(prompt) && xyradius_vol_cut.isIn(prompt))
        ) {
            LogInfo << "Prompt is a chimney\n";
            continue;
        }

        if (prompt.totq < prompt_lower_thold || prompt_upper_thold < prompt.totq) {
            LogInfo << "Prompt not in energy range\n";
            continue;
        }

        bool is_vetoed = false;
        for (const WaterPoolMuonVetoSelection& cut : mu_wp_bundle_cut) {
            if (!cut.isIn(prompt)) continue;
            is_vetoed = true;
            break;
        }
        if (is_vetoed) {
            LogInfo << "Prompt is muon vetoed\n";
            continue;
        }

        is_vetoed = false;
        const track* trk_cosmo = nullptr;
        for (const BasicMuonVetoSelection& cut : mu_cosmo_cut) {
            if (!cut.isIn(prompt)) continue;
            is_vetoed = true;
            trk_cosmo = &cut.m_trk;
            break;
        }
        if (!is_vetoed) {
            LogInfo << "Prompt is not in cylindrical muon cosmogenic cut\n";
            continue;
        }

        WindowTimeSelection multi_prompt_time{prompt.ts, TimeStamp{0, -2000000}, TimeStamp{0, 0}};
        bool prompt_has_multi = false;
        for (const vertex& cand : bef_vertices) {
            if (!multi_prompt_time.isIn(cand)) continue;
            if (!fiducial_vol_cut.isIn(cand)) continue;
            if (
                (upper_height_vol_cut.isIn(cand) && xyradius_vol_cut.isIn(cand)) ||
                (lower_height_vol_cut.isIn(cand) && xyradius_vol_cut.isIn(cand))
            ) continue;
            if (cand.totq < prompt_lower_thold || prompt_upper_thold < cand.totq) continue;
            
            is_vetoed = false;
            for (const WaterPoolMuonVetoSelection& cut : mu_wp_bundle_cut) {
                if (!cut.isIn(cand)) continue;
                is_vetoed = true;
                break;
            }
            if (is_vetoed) continue;
            is_vetoed = false;
            for (const BasicMuonVetoSelection& cut : mu_cosmo_cut) {
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

        WindowTimeSelection correlation_time_cut{prompt.ts, TimeStamp{0, 5000}, TimeStamp{0, 2000000}};
        SphereVolumeSelection distance_correlation_cut{prompt.pos, 1500.0};

        for (const vertex& delayed : aft_vertices) {
            LogInfo << delayed << '\n';
            if (!fiducial_vol_cut.isIn(delayed)) {
                LogInfo << "Delayed not in fiducial volume\n";
                continue;
            }

            if (
                (upper_height_vol_cut.isIn(delayed) && xyradius_vol_cut.isIn(delayed)) ||
                (lower_height_vol_cut.isIn(delayed) && xyradius_vol_cut.isIn(delayed))
            ) {
                LogInfo << "Delayed is a chimney\n";
                continue;
            }

            if (delayed.totq < delayed_lower_thold || delayed_upper_thold < delayed.totq) {
                LogInfo << "Delayed not in energy range\n";
                continue;
            }

            if (!correlation_time_cut.isIn(delayed)) {
                LogInfo << "Delayed is not correlated in time\n";
                continue;
            }
            if (!distance_correlation_cut.isIn(delayed)) {
                LogInfo << "Delayed is not correlated in space\n";
                continue;
            }

            is_vetoed = false;
            for (const WaterPoolMuonVetoSelection& cut : mu_wp_bundle_cut) {
                if (!cut.isIn(delayed)) continue;
                is_vetoed = true;
                break;
            }
            if (is_vetoed) {
                LogInfo << "Delayed is muon vetoed\n";
                continue;
            }

            is_vetoed = false;
            for (const BasicMuonVetoSelection& cut : mu_cosmo_cut) {
                if (!cut.isIn(delayed)) continue;
                is_vetoed = true;
                break;
            }
            if (!is_vetoed) {
                LogInfo << "Delayed is not in cylindrical muon cosmogenic cut\n";
                continue;
            }

            WindowTimeSelection multi_delayed_time{delayed.ts, TimeStamp{0, 0}, TimeStamp{0, 2000000}};
            bool delayed_has_multi = false;
            for (const vertex& cand : aft_vertices) {
                if (cand.ts == delayed.ts) continue; // same event
                if (!fiducial_vol_cut.isIn(cand)) continue;
                if (
                    (upper_height_vol_cut.isIn(cand) && xyradius_vol_cut.isIn(cand)) ||
                    (lower_height_vol_cut.isIn(cand) && xyradius_vol_cut.isIn(cand))
                ) continue;
                if (cand.totq < prompt_lower_thold || prompt_upper_thold < cand.totq) continue;
                
                is_vetoed = false;
                for (const WaterPoolMuonVetoSelection& cut : mu_wp_bundle_cut) {
                    if (!cut.isIn(cand)) continue;
                    is_vetoed = true;
                    break;
                }
                if (is_vetoed) continue;
                is_vetoed = false;
                for (const BasicMuonVetoSelection& cut : mu_cosmo_cut) {
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
            tracks_for_cosmo.push_back(*trk_cosmo);
            LogInfo << "IBD event detected!\n";
        }
    }

    for (std::size_t k = 0ul; k < cosmos.size(); ++k) {
        const ibd& ibd_ = cosmos[k];
        const track& trk_ = tracks_for_cosmo[k];
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