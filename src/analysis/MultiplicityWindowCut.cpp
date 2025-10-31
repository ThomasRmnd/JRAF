#include "analysis/MultiplicityWindowCut.hpp"

#include <algorithm>

#include "SniperKernel/SniperLog.h"

#include "event/Event.hpp"
#include "event/EventCache.hpp"
#include "event/IBD.hpp"
#include "selection/Energy.hpp"
#include "selection/Muon.hpp"
#include "selection/Vertex.hpp"
#include "selection/Volume.hpp"

MultiplicityWindowCut::MultiplicityWindowCut(const std::string& name, const std::string& method) : 
    Analysis{name, method} 
{}

bool MultiplicityWindowCut::initialize() {
    if (!Analysis::initialize()) return false;
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
    extractEvent(buf, tracks, cur_vertices, bef_vertices, aft_vertices);


    std::vector<TimeRangeMuonVetoSelection> mu_cut;
    for (std::vector<std::vector<track>>::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
        if (it->empty()) continue;
        mu_cut.emplace_back(it->front(), TimeStamp{0, 0}, TimeStamp{0, 5000000});
    }

    FiducialVolumeSelection fiducial_vol_cut{16500.0};
    ChimneySelection chimney_cut{15500.0, 3000.0};
    // ChargeRangeSelection prompt_charge_cut{1500.0, 20000.0};
    // ChargeRangeSelection delayed_charge_cut{4000.0, 6000.0};
    EnergyRangeSelection prompt_energy_cut{0.7, 12.0};
    EnergyRangeSelection delayed_energy_cut{2.0, 2.5};
    EnergyRangeSelection multiplicity_energy_cut{2.0, 12.0};

    std::vector<ibd> ibds;
    std::vector<vertex> multi_vertices;
    std::vector<unsigned char> multi_window_types;

    for (const vertex& prompt : cur_vertices) {
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
        for (const TimeRangeMuonVetoSelection& cut : mu_cut) {
            if (!cut.isIn(prompt)) continue;
            is_vetoed = true;
            break;
        }
        if (is_vetoed) {
            LogInfo << "Prompt is muon vetoed\n";
            continue;
        }

        TimeRangeSelection multi_prompt_time{prompt.ts, TimeStamp{0, -1000000}, TimeStamp{0, 0}};
        bool prompt_has_multi = false;
        for (const vertex& cand : bef_vertices) {
            if (!multi_prompt_time.isIn(cand)) continue;
            // if (!fiducial_vol_cut.isIn(cand)) continue;
            // if (chimney_cut.isIn(cand)) continue;
            if (!multiplicity_energy_cut.isIn(cand)) continue;
            is_vetoed = false;
            for (const TimeRangeMuonVetoSelection& cut : mu_cut) {
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

        VertexCorrelationSelection correlation_cut{prompt, 1500.0, TimeStamp{0, 5000}, TimeStamp{0, 1000000}};

        for (const vertex& delayed : aft_vertices) {
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
            for (const TimeRangeMuonVetoSelection& cut : mu_cut) {
                if (!cut.isIn(delayed)) continue;
                is_vetoed = true;
                break;
            }
            if (is_vetoed) {
                LogInfo << "Delayed is muon vetoed\n";
                continue;
            }

            TimeRangeSelection multi_delayed_time{delayed.ts, TimeStamp{0, 0}, TimeStamp{0, 1000000}};
            bool delayed_has_multi = false;
            for (const vertex& cand : aft_vertices) {
                if (cand.ts == delayed.ts) continue; // same event
                // if (!fiducial_vol_cut.isIn(cand)) continue;
                // if (chimney_cut.isIn(cand)) continue;
                if (!multiplicity_energy_cut.isIn(cand)) continue;
                is_vetoed = false;
                for (const TimeRangeMuonVetoSelection& cut : mu_cut) {
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