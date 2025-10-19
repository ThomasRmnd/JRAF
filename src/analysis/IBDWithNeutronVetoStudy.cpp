#include "analysis/IBDWithNeutronVetoStudy.hpp"

#include <algorithm>
#include <chrono>

#include "SniperKernel/SniperLog.h"

#include "event/Event.hpp"
#include "event/EventCache.hpp"
#include "event/IBD.hpp"
#include "selection/Energy.hpp"
#include "selection/Muon.hpp"
#include "selection/Volume.hpp"

IBDWithNeutronVetoStudy::IBDWithNeutronVetoStudy(const std::string& name) : 
    Analysis{name} 
{}

bool IBDWithNeutronVetoStudy::initialize() {
    if (!Analysis::initialize()) return false;
    m_tree->Branch("totq_p", &totq_p);
    m_tree->Branch("totq_d", &totq_d);
    return true; 
}

void IBDWithNeutronVetoStudy::process(JM::NavBuffer* buf) {
    std::vector<std::vector<track>> tracks;
    std::vector<vertex> cur_vertices;
    std::vector<vertex> bef_vertices;
    std::vector<vertex> aft_vertices;
    tracks.reserve(buf->size());
    bef_vertices.reserve(buf->size() / 2);
    aft_vertices.reserve(buf->size() / 2);

    for (JM::NavBuffer::Iterator it = buf->begin(); it != buf->end(); ++it) {
        JM::EvtNavigator* nav = it->get();
        if (!nav) continue;

        std::shared_ptr<Event> evt_ptr = EventCache::load(nav);
        if (!evt_ptr) continue;

        const Event& evt = *evt_ptr;

        tracks.push_back(evt.tracks);
        if (it < buf->current()) {
            bef_vertices.insert(bef_vertices.end(), evt.vertices.begin(), evt.vertices.end());
        } else if (buf->current() < it) {
            aft_vertices.insert(aft_vertices.end(), evt.vertices.begin(), evt.vertices.end());
        } else {
            cur_vertices.insert(cur_vertices.end(), evt.vertices.begin(), evt.vertices.end());
        }
    }

    std::vector<WaterPoolMuonVetoSelection> mu_cut;
    for (std::vector<std::vector<track>>::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
        if (it->empty()) continue;
        mu_cut.emplace_back(it->front(), TimeStamp{0, 5000000});
    }

    FiducialVolumeSelection fiducial_vol_cut{17200.0};
    HeightVolumeSelection lower_height_vol_cut{-20050.0, -11000.0};
    HeightVolumeSelection upper_height_vol_cut{ 11000.0,  20050.0};
    XYRadiusVolumeSelection xyradius_vol_cut{0.0, 3000.0};
    double prompt_lower_thold = 1500.0;
    double prompt_upper_thold = 20000.0;
    double delayed_lower_thold = 3700.0;
    double delayed_upper_thold = 6000.0;

    double spa_neu_radius = 3000.0;
    TimeStamp spa_neu_time_window = TimeStamp{0, 1500000000};
    std::vector<SphereVolumeSelection> spa_neu_cut_sph;
    std::vector<UpperTimeSelection> spa_neu_cut_time;

    for (const vertex& neu : bef_vertices) {
        bool is_in_veto = false;
        for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
            if (!cut.isIn(neu)) continue;
            is_in_veto = true;
            break;
        }
        if (!is_in_veto) continue;
        if (neu.totq < prompt_lower_thold || prompt_upper_thold < neu.totq) continue;
        spa_neu_cut_sph.emplace_back(neu.pos, spa_neu_radius);
        spa_neu_cut_time.emplace_back(neu.ts, spa_neu_time_window);
    }
    for (const vertex& neu : cur_vertices) {
        bool is_in_veto = false;
        for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
            if (!cut.isIn(neu)) continue;
            is_in_veto = true;
            break;
        }
        if (!is_in_veto) continue;
        if (neu.totq < prompt_lower_thold || prompt_upper_thold < neu.totq) continue;
        spa_neu_cut_sph.emplace_back(neu.pos, spa_neu_radius);
        spa_neu_cut_time.emplace_back(neu.ts, spa_neu_time_window);
    }
    for (const vertex& neu : aft_vertices) {
        bool is_in_veto = false;
        for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
            if (!cut.isIn(neu)) continue;
            is_in_veto = true;
            break;
        }
        if (!is_in_veto) continue;
        if (neu.totq < prompt_lower_thold || prompt_upper_thold < neu.totq) continue;
        spa_neu_cut_sph.emplace_back(neu.pos, spa_neu_radius);
        spa_neu_cut_time.emplace_back(neu.ts, spa_neu_time_window);
    }

    std::vector<ibd> ibds;

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
        for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
            if (!cut.isIn(prompt)) continue;
            is_vetoed = true;
            break;
        }
        if (is_vetoed) {
            LogInfo << "Prompt is muon vetoed\n";
            continue;
        }

        bool is_prompt_neu_spa = false;
        for (std::size_t k = 0ul; k < spa_neu_cut_sph.size(); ++k) {
            if (!spa_neu_cut_sph[k].isIn(prompt) || !spa_neu_cut_time[k].isIn(prompt)) continue;
            is_prompt_neu_spa = true;
        }
        if (is_prompt_neu_spa) {
            LogInfo << "Prompt is in spatial-temporal neutron veto\n";
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
            for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
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
                LogInfo << "Delayed not correlated in time\n";
                continue;
            }
            if (!distance_correlation_cut.isIn(delayed)) {
                LogInfo << "Delayed not correlated in space\n";
                continue;
            }

            is_vetoed = false;
            for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
                if (!cut.isIn(delayed)) continue;
                is_vetoed = true;
                break;
            }
            if (is_vetoed) {
                LogInfo << "Delayed is muon vetoed\n";
                continue;
            }

            bool is_delayed_neu_spa = false;
            for (std::size_t k = 0ul; k < spa_neu_cut_sph.size(); ++k) {
                if (!spa_neu_cut_sph[k].isIn(delayed) || !spa_neu_cut_time[k].isIn(delayed)) continue;
                is_delayed_neu_spa = true;
            }
            if (is_delayed_neu_spa) {
                LogInfo << "Delayed is in spatial-temporal neutron veto\n";
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
                for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
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

            ibds.emplace_back(prompt, delayed);
            LogInfo << "IBD event detected!\n";
        }
    }

    for (const ibd& ibd_ : ibds) {
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
        m_tree->Fill();
    }
}