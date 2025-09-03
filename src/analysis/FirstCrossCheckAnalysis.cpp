#include "analysis/FirstCrossCheckAnalysis.hpp"

#include <algorithm>

#include "SniperKernel/SniperLog.h"

#include "event/Event.hpp"
#include "event/IBD.hpp"
#include "selection/Energy.hpp"
#include "selection/Muon.hpp"
#include "selection/Volume.hpp"

FirstCrossCheckAnalysis::FirstCrossCheckAnalysis(const std::string& name) : 
    Analysis{name} 
{}

bool FirstCrossCheckAnalysis::initialize() {
    return true;
}

void FirstCrossCheckAnalysis::process(JM::NavBuffer* buf) {
    Event __evt;
    __evt.load(buf->curEvt());
    LogInfo << __evt << '\n';

    std::vector<std::vector<track>> tracks;
    std::vector<vertex> cur_vertices;
    std::vector<vertex> bef_vertices;
    std::vector<vertex> aft_vertices;
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
            cur_vertices.insert(cur_vertices.end(), evt.vertices.begin(), evt.vertices.end());
        }
    }

    std::vector<WaterPoolMuonVetoSelection> mu_cut;
    for (std::vector<std::vector<track>>::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
        if (it->empty()) continue;
        mu_cut.emplace_back(it->front(), TimeStamp{0, 2000000});
    }

    FiducialVolumeSelection fiducial_vol_cut{17200.0};
    HeightVolumeSelection height_vol_cut{-17700.0, 11000.0};
    XYRadiusVolumeSelection xyradius_vol_cut{3000.0, 17700.0};
    double prompt_lower_thold = 1500.0;
    double prompt_upper_thold = 20000.0;
    double delayed_lower_thold = 3700.0;
    double delayed_upper_thold = 6000.0;

    std::vector<ibd> ibds;

    for (const vertex& prompt : cur_vertices) {
        LogInfo << prompt << '\n';
        if (!fiducial_vol_cut.isIn(prompt)) continue;
        LogInfo << "Prompt in fiducial volume\n";
        if (!height_vol_cut.isIn(prompt)) continue;
        LogInfo << "Prompt below max height\n";
        if (!xyradius_vol_cut.isIn(prompt)) continue;
        LogInfo << "Prompt outside XY radius\n";

        if (prompt.totq < prompt_lower_thold || prompt_upper_thold < prompt.totq) continue;
        LogInfo << "Prompt in energy range\n";

        bool is_vetoed = false;
        for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
            LogInfo << prompt.ts << " - " << cut.trkptr->ts << " = " << prompt.ts - cut.trkptr->ts << ' ' << (prompt.ts - cut.trkptr->ts < TimeStamp{0, 2000000}) << ' ' << (!cut.isIn(prompt)) << '\n';
            if (!cut.isIn(prompt)) continue;
            is_vetoed = true;
            break;
        }
        if (is_vetoed) continue;
        LogInfo << "Prompt is not vetoed\n";

        WindowTimeSelection multi_prompt_time{prompt.ts, TimeStamp{0, -2000000}, TimeStamp{0, 0}};
        bool prompt_has_multi = false;
        for (const vertex& cand : bef_vertices) {
            if (!multi_prompt_time.isIn(cand)) continue;
            prompt_has_multi = true;
            break;
        }
        if (prompt_has_multi) continue;
        LogInfo << "Prompt has no multiplicity\n";

        WindowTimeSelection correlation_time_cut{prompt.ts, TimeStamp{0, 1000}, TimeStamp{0, 2000000}};
        SphereVolumeSelection distance_correlation_cut{prompt.pos, 1500.0};

        for (const vertex& delayed : aft_vertices) {
            if (!fiducial_vol_cut.isIn(delayed)) continue;
            if (!height_vol_cut.isIn(delayed)) continue;
            if (!xyradius_vol_cut.isIn(delayed)) continue;

            if (delayed.totq < delayed_lower_thold || delayed_upper_thold < delayed.totq) continue;

            if (!correlation_time_cut.isIn(delayed)) continue;
            if (!distance_correlation_cut.isIn(delayed)) continue;

            is_vetoed = false;
            for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
                if (!cut.isIn(delayed)) continue;
                is_vetoed = true;
                break;
            }
            if (is_vetoed) continue;

            WindowTimeSelection multi_delayed_time{delayed.ts, TimeStamp{0, 0}, TimeStamp{0, 2000000}};
            bool delayed_has_multi = false;
            for (const vertex& cand : aft_vertices) {
                if (cand.ts <= delayed.ts) continue;
                if (!multi_delayed_time.isIn(cand)) continue;
                delayed_has_multi = true;
                break;
            }
            if (delayed_has_multi) continue;

            ibds.emplace_back(prompt, delayed);
        }
    }

    for (const ibd& ibd_ : ibds) {
        posx_p = ibd_.prompt.pos.x;
        posy_p = ibd_.prompt.pos.y;
        posz_p = ibd_.prompt.pos.z;
        e_p = ibd_.prompt.energy;
        sec_p = ibd_.prompt.ts.GetSec();
        nsec_p = ibd_.prompt.ts.GetNanoSec();
        posx_d = ibd_.delayed.pos.x;
        posy_d = ibd_.delayed.pos.y;
        posz_d = ibd_.delayed.pos.z;
        e_d = ibd_.delayed.energy;
        sec_d = ibd_.delayed.ts.GetSec();
        nsec_d = ibd_.delayed.ts.GetNanoSec();
        m_tree->Fill();
    }
}