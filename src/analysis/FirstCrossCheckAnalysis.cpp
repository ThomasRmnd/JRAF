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
    if (!Analysis::initialize()) return false;
    m_tree->Branch("totq_p", &totq_p);
    m_tree->Branch("totq_d", &totq_d);
    m_tree->Branch("daq_sec", &m_daq_sec);
    m_tree->Branch("daq_nsec", &m_daq_nsec);
    m_tree->Branch("muveto_sec", &m_muveto_sec);
    m_tree->Branch("muveto_nsec", &m_muveto_nsec);
    return true; 
}

void FirstCrossCheckAnalysis::process(JM::NavBuffer* buf) {
    Event __evt;
    __evt.load(buf->curEvt());
    LogInfo << __evt << '\n';

    TimeStamp daq_ts{m_daq_sec, m_daq_nsec};
    daq_ts.Add(__evt.ts - m_prv_ts);
    m_daq_sec = daq_ts.GetSec();
    m_daq_nsec = daq_ts.GetNanoSec();

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
        if (it == buf->current() && !evt.tracks.empty()) {
            TimeStamp muveto_ts{m_muveto_sec, m_muveto_nsec};
            muveto_ts.Add(TimeStamp{0, 2000000});
            m_muveto_sec = muveto_ts.GetSec();
            m_muveto_nsec = muveto_ts.GetNanoSec();
        }
    }


    std::vector<WaterPoolMuonVetoSelection> mu_cut;
    for (std::vector<std::vector<track>>::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
        if (it->empty()) continue;
        mu_cut.emplace_back(it->front(), TimeStamp{0, 2000000});
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

    for (const vertex& prompt : cur_vertices) {
        LogInfo << prompt << '\n';
        if (!fiducial_vol_cut.isIn(prompt)) continue;
        LogInfo << "Prompt in fiducial volume\n";
        if (
            (upper_height_vol_cut.isIn(prompt) && xyradius_vol_cut.isIn(prompt)) ||
            (lower_height_vol_cut.isIn(prompt) && xyradius_vol_cut.isIn(prompt))
        ) continue;
        LogInfo << "Prompt is not a chimney\n";

        if (prompt.totq < prompt_lower_thold || prompt_upper_thold < prompt.totq) continue;
        LogInfo << "Prompt in energy range\n";

        bool is_vetoed = false;
        for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
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
            if (!fiducial_vol_cut.isIn(cand)) continue;
            is_vetoed = false;
            for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
                if (!cut.isIn(cand)) continue;
                is_vetoed = true;
                break;
            }
            if (is_vetoed) continue;
            if (cand.totq < prompt_lower_thold || prompt_upper_thold < cand.totq) continue;
            prompt_has_multi = true;
            break;
        }
        if (prompt_has_multi) continue;
        LogInfo << "Prompt has no multiplicity\n";

        WindowTimeSelection correlation_time_cut{prompt.ts, TimeStamp{0, 5000}, TimeStamp{0, 2000000}};
        SphereVolumeSelection distance_correlation_cut{prompt.pos, 1500.0};

        for (const vertex& delayed : aft_vertices) {
            LogInfo << delayed << '\n';
            if (!fiducial_vol_cut.isIn(delayed)) continue;
            LogInfo << "Delayed in fiducial volume\n";
            if (
                (upper_height_vol_cut.isIn(prompt) && xyradius_vol_cut.isIn(prompt)) ||
                (lower_height_vol_cut.isIn(prompt) && xyradius_vol_cut.isIn(prompt))
            ) continue;
            LogInfo << "Delayed is not a chimney\n";

            if (delayed.totq < delayed_lower_thold || delayed_upper_thold < delayed.totq) continue;
            LogInfo << "Delayed in energy range\n";

            if (!correlation_time_cut.isIn(delayed)) continue;
            LogInfo << "Delayed is correlated in time\n";
            if (!distance_correlation_cut.isIn(delayed)) continue;
            LogInfo << "Delayed is correlated in space\n";

            is_vetoed = false;
            for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
                if (!cut.isIn(delayed)) continue;
                is_vetoed = true;
                break;
            }
            if (is_vetoed) continue;
            LogInfo << "Delayed is is not vetoed\n";

            WindowTimeSelection multi_delayed_time{delayed.ts, TimeStamp{0, 0}, TimeStamp{0, 2000000}};
            bool delayed_has_multi = false;
            for (const vertex& cand : aft_vertices) {
                if (cand.ts == delayed.ts) continue; // same event
                if (!fiducial_vol_cut.isIn(cand)) continue;
                is_vetoed = false;
                for (const WaterPoolMuonVetoSelection& cut : mu_cut) {
                    if (!cut.isIn(cand)) continue;
                    is_vetoed = true;
                    break;
                }
                if (is_vetoed) continue;
                if (cand.ts < delayed.ts && (prompt_lower_thold <= cand.totq && cand.totq <= prompt_upper_thold)) {
                    delayed_has_multi = true;
                    break; // in-between p-d multiplicity
                }
                if (!multi_delayed_time.isIn(cand)) continue;
                if (cand.totq < prompt_lower_thold || prompt_upper_thold < cand.totq) continue;
                delayed_has_multi = true;
                break; // after p-d multiplicity
            }
            if (delayed_has_multi) continue;
            LogInfo << "Delayed has no multiplicity\n";

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