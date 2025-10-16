#include "analysis/MultiplicityWindowCut.hpp"

#include <algorithm>

#include "SniperKernel/SniperLog.h"

#include "event/Event.hpp"
#include "event/IBD.hpp"
#include "selection/Energy.hpp"
#include "selection/Muon.hpp"
#include "selection/Volume.hpp"

MultiplicityWindowCut::MultiplicityWindowCut(const std::string& name) : 
    Analysis{name} 
{}

bool MultiplicityWindowCut::initialize() {
    if (!Analysis::initialize()) return false;
    m_tree->Branch("totq_p", &totq_p);
    m_tree->Branch("totq_d", &totq_d);

    m_tree->Branch("window_type", &m_window_type);
    m_tree->Branch("posx_m", &posx_m);
    m_tree->Branch("posy_m", &posy_m);
    m_tree->Branch("posz_m", &posz_m);
    m_tree->Branch("e_m", &e_m);
    m_tree->Branch("totq_m", &totq_m);
    m_tree->Branch("sec_m", &sec_m);
    m_tree->Branch("nsec_m", &nsec_m);
    return true; 
}

void MultiplicityWindowCut::process(JM::NavBuffer* buf) {
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

    std::vector<ibd> ibds;
    std::vector<vertex> multi_vertices;
    std::vector<unsigned char> multi_window_types;

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
            multi_vertices.push_back(cand);
            multi_window_types.push_back(0);
            LogInfo << "An event is in the prompt multiplicity cut by " << cand.ts << '\n';
            break;
        }
        // if (prompt_has_multi) continue;

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
                    multi_vertices.push_back(cand);
                    multi_window_types.push_back(1);
                    LogInfo << "An event is in the in-between multiplicity cut by " << cand.ts << '\n';
                    break; // in-between p-d multiplicity
                }
                if (!multi_delayed_time.isIn(cand)) continue;
                delayed_has_multi = true;
                multi_vertices.push_back(cand);
                multi_window_types.push_back(2);
                LogInfo << "An event is in the after multiplicity cut by " << cand.ts << '\n';
                break; // after p-d multiplicity
            }
            // if (delayed_has_multi) continue;

            if (!prompt_has_multi || !delayed_has_multi) continue;
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
        for (std::size_t k = 0; k < multi_vertices.size(); ++k) {
            const vertex& v = multi_vertices[k];
            m_window_type = multi_window_types[k];
            posx_m = v.pos.x;
            posy_m = v.pos.y;
            posz_m = v.pos.z;
            e_m = v.energy;
            totq_m = v.totq;
            sec_m = v.ts.GetSec();
            nsec_m = v.ts.GetNanoSec();
            m_tree->Fill();
        }
    }
}