#include "analysis/IBDWithNeutronVetoStudy.hpp"

#include <algorithm>
#include <chrono>

#include "SniperKernel/SniperLog.h"

#include "analysis/NavBufferCache.hpp"
#include "event/Event.hpp"
#include "event/EventCache.hpp"
#include "event/IBD.hpp"
#include "selection/Energy.hpp"
#include "selection/Muon.hpp"
#include "selection/Vertex.hpp"
#include "selection/Volume.hpp"

IBDWithNeutronVetoStudy::IBDWithNeutronVetoStudy(const std::string& name, const std::string& method) : 
    Analysis{name, method} 
{}

void IBDWithNeutronVetoStudy::process(JM::NavBuffer*) {
    const std::vector<std::vector<track>>& tracks = NavBufferCache::getTracks(m_method);
    const std::vector<vertex>& cur_vertices = NavBufferCache::getVertices(m_method, NavBufferCache::VertexRegion::Current);
    const std::vector<vertex>& bef_vertices = NavBufferCache::getVertices(m_method, NavBufferCache::VertexRegion::Before);
    const std::vector<vertex>& aft_vertices = NavBufferCache::getVertices(m_method, NavBufferCache::VertexRegion::After);
    // extractEvent(buf, tracks, cur_vertices, bef_vertices, aft_vertices);

    std::vector<TimeRangeMuonVetoSelection> mu_cut;
    std::vector<TimeRangeMuonVetoSelection> mu_spa_neu_cut;
    for (std::vector<std::vector<track>>::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
        if (it->empty()) continue;
        mu_cut.emplace_back(it->front(), TimeStamp{0, 0}, TimeStamp{0, 5000000});
        mu_spa_neu_cut.emplace_back(it->front(), TimeStamp{0, 20000}, TimeStamp{0, 2000000});
    }

    FiducialVolumeSelection fiducial_vol_cut{16500.0};
    ChimneySelection chimney_cut{15500.0, 3000.0};
    // ChargeRangeSelection prompt_charge_cut{1500.0, 20000.0};
    // ChargeRangeSelection delayed_charge_cut{4000.0, 6000.0};
    EnergyRangeSelection prompt_energy_cut{0.7, 12.0};
    EnergyRangeSelection delayed_energy_cut{2.0, 2.5};
    EnergyRangeSelection spa_neu_energy_cut{1.5, 20.0};
    EnergyRangeSelection multiplicity_energy_cut{2.0, 12.0};

    std::vector<VertexCorrelationSelection> spa_neu_cut;
    for (const vertex& neu : bef_vertices) {
        bool is_in_veto = false;
        for (const TimeRangeMuonVetoSelection& cut : mu_spa_neu_cut) {
            if (!cut.isIn(neu)) continue;
            is_in_veto = true;
            break;
        }
        if (!is_in_veto) continue;
        if (!spa_neu_energy_cut.isIn(neu)) continue;
        spa_neu_cut.emplace_back(neu, 4000.0, TimeStamp{0, 0}, TimeStamp{0, 1200000000});
    }
    for (const vertex& neu : cur_vertices) {
        bool is_in_veto = false;
        for (const TimeRangeMuonVetoSelection& cut : mu_spa_neu_cut) {
            if (!cut.isIn(neu)) continue;
            is_in_veto = true;
            break;
        }
        if (!is_in_veto) continue;
        if (!spa_neu_energy_cut.isIn(neu)) continue;
        spa_neu_cut.emplace_back(neu, 4000.0, TimeStamp{0, 0}, TimeStamp{0, 1200000000});
    }
    for (const vertex& neu : aft_vertices) {
        bool is_in_veto = false;
        for (const TimeRangeMuonVetoSelection& cut : mu_spa_neu_cut) {
            if (!cut.isIn(neu)) continue;
            is_in_veto = true;
            break;
        }
        if (!is_in_veto) continue;
        if (!spa_neu_energy_cut.isIn(neu)) continue;
        spa_neu_cut.emplace_back(neu, 4000.0, TimeStamp{0, 0}, TimeStamp{0, 1200000000});
    }

    std::vector<ibd> ibds;

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

        bool is_prompt_neu_spa = false;
        for (const VertexCorrelationSelection& cut : spa_neu_cut) {
            if (!cut.isIn(prompt)) continue;
            is_prompt_neu_spa = true;
            break;
        }
        if (is_prompt_neu_spa) {
            LogInfo << "Prompt is in spatial-temporal neutron veto\n";
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
            break;
        }
        if (prompt_has_multi) {
            LogInfo << "Prompt has multiplicity\n";
            continue;
        }

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

            bool is_delayed_neu_spa = false;
            for (const VertexCorrelationSelection& cut : spa_neu_cut) {
                if (!cut.isIn(delayed)) continue;
                is_delayed_neu_spa = true;
                break;
            }
            if (is_delayed_neu_spa) {
                LogInfo << "Delayed is in spatial-temporal neutron veto\n";
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