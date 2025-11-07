#include "analysis/CosmoRateWithNeutronAnalysis.hpp"

#include <set>

#include "SniperKernel/SniperLog.h"

#include "event/IBD.hpp"
#include "selection/Energy.hpp"
#include "selection/Vertex.hpp"
#include "selection/Volume.hpp"

CosmoRateWithNeutronAnalysis::CosmoRateWithNeutronAnalysis(const std::string& name, const std::string& method) :
    Analysis{name, method}
{}

bool CosmoRateWithNeutronAnalysis::initialize() {
    if (!Analysis::initialize()) return false;
    m_tree->Branch("dlat_p", &dlat_p); 
    m_tree->Branch("dlat_d", &dlat_d);
    m_tree->Branch("dt_mu2p_sec", &dt_mu2p_sec); 
    m_tree->Branch("dt_mu2d_sec", &dt_mu2d_sec);
    m_tree->Branch("dt_mu2p_nsec", &dt_mu2p_nsec), 
    m_tree->Branch("dt_mu2d_nsec", &dt_mu2d_nsec);
    m_tree->Branch("posx_n", &posx_n);
    m_tree->Branch("posy_n", &posy_n);
    m_tree->Branch("posz_n", &posz_n);
    m_tree->Branch("e_n", &e_n);
    m_tree->Branch("sec_n", &sec_n);
    m_tree->Branch("nsec_n", &nsec_n);
    m_tree->Branch("totq_n", &totq_n);
    return true;
}

void CosmoRateWithNeutronAnalysis::process(const EventContext::View& events) {
    std::vector<TimeRangeMuonVetoSelection> mu_cut;
    std::vector<MuonAssociatedWithNeutron> mu_neu_cut;
    std::set<TimeStamp> visited;
    for (const track& trk : events.tracks()) {
        if (visited.find(trk.ts) != visited.end()) continue;
        visited.insert(trk.ts);
        mu_cut.emplace_back(trk, TimeStamp{0, 0}, TimeStamp{0, 5000000});
        mu_neu_cut.emplace_back(
            TimeRangeMuonVetoSelection{trk, TimeStamp{0, 5000000}, TimeStamp{0, 1200000000}},
            TimeRangeMuonVetoSelection{trk, TimeStamp{0, 20000}, TimeStamp{0, 2000000}}
        );
    }

    FiducialVolumeSelection fiducial_vol_cut{16500.0};
    ChimneySelection chimney_cut{15500.0, 3000.0};
    // ChargeRangeSelection prompt_charge_cut{1500.0, 20000.0};
    // ChargeRangeSelection delayed_charge_cut{4000.0, 6000.0};
    EnergyRangeSelection prompt_energy_cut{0.7, 12.0};
    EnergyRangeSelection delayed_energy_cut{2.0, 2.5};
    EnergyRangeSelection multiplicity_energy_cut{2.0, 12.0};

    for (const vertex& neu : events.vertices()) {
        if (!delayed_energy_cut.isIn(neu)) continue;
        for (MuonAssociatedWithNeutron& mu : mu_neu_cut) {
            if (!mu.neu_veto.isIn(neu)) continue;
            mu.neu.push_back(neu);
        }
    }

    std::vector<ibd> ibds;

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
        for (const vertex& cand : events.before()) {
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
            for (const vertex& cand : events.after()) {
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
        for (const MuonAssociatedWithNeutron& mu : mu_neu_cut) {
            if (!mu.cosmo_veto.isIn(ibd_.prompt) || !mu.cosmo_veto.isIn(ibd_.delayed)) continue;
            if (mu.neu.empty()) {
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

                const track& trk_ = mu.cosmo_veto.c_trk;
                vec3 trk_dir = unit(trk_.fpos - trk_.ipos);
                dlat_p = mag(cross(trk_dir, ibd_.prompt.pos - trk_.ipos));
                dlat_d = mag(cross(trk_dir, ibd_.delayed.pos - trk_.ipos));
                dt_mu2p_sec = (ibd_.prompt.ts - trk_.ts).GetSec();
                dt_mu2d_sec = (ibd_.delayed.ts - trk_.ts).GetSec();
                dt_mu2p_nsec = (ibd_.prompt.ts - trk_.ts).GetNanoSec();
                dt_mu2d_nsec = (ibd_.delayed.ts - trk_.ts).GetNanoSec();

                posx_n = 0.0;
                posy_n = 0.0;
                posz_n = 0.0;
                e_n = 0.0;
                sec_n = 0l;
                nsec_n = 0;
                totq_n = 0.0;

                m_tree->Fill();
            }
            else {
                for (const vertex& neu : mu.neu) {
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

                    const track& trk_ = mu.cosmo_veto.c_trk;
                    vec3 trk_dir = unit(trk_.fpos - trk_.ipos);
                    dlat_p = mag(cross(trk_dir, ibd_.prompt.pos - trk_.ipos));
                    dlat_d = mag(cross(trk_dir, ibd_.delayed.pos - trk_.ipos));
                    dt_mu2p_sec = (ibd_.prompt.ts - trk_.ts).GetSec();
                    dt_mu2d_sec = (ibd_.delayed.ts - trk_.ts).GetSec();
                    dt_mu2p_nsec = (ibd_.prompt.ts - trk_.ts).GetNanoSec();
                    dt_mu2d_nsec = (ibd_.delayed.ts - trk_.ts).GetNanoSec();

                    posx_n = neu.pos.x;
                    posy_n = neu.pos.y;
                    posz_n = neu.pos.z;
                    e_n = neu.energy;
                    sec_n = neu.ts.GetSec();
                    nsec_n = neu.ts.GetNanoSec();
                    totq_n = neu.totq;

                    m_tree->Fill();
                }
            }
        }
        
    }
}