#include "analysis/IBDAnalysis.hpp"

#include <algorithm>

#include "SniperKernel/SniperLog.h"

#include "event/IBD.hpp"
#include "selection/Energy.hpp"
#include "selection/Muon.hpp"
#include "selection/Vertex.hpp"
#include "selection/Volume.hpp"

IBDAnalysis::IBDAnalysis(const std::string& name, const std::string& method) : 
    Analysis{name, method} 
{}

bool IBDAnalysis::initialize() {
    if (!Analysis::initialize()) return false;
    
    m_tree->Branch("posx_n", &posx_n);
    m_tree->Branch("posy_n", &posy_n);
    m_tree->Branch("posz_n", &posz_n);
    m_tree->Branch("e_n", &e_n);
    m_tree->Branch("sec_n", &sec_n);
    m_tree->Branch("nsec_n", &nsec_n);

    m_tree->Branch("totq_n", &totq_n);
    m_tree->Branch("meanq_n", &meanq_n);
    m_tree->Branch("stdq_n", &stdq_n);
    m_tree->Branch("minq_n", &minq_n);
    m_tree->Branch("maxq_n", &maxq_n);
    m_tree->Branch("nhit_n", &nhit_n);
    m_tree->Branch("meant_n", &meant_n);
    m_tree->Branch("stdt_n", &stdt_n);

    m_tree->Branch("posx_mult", &posx_mult);
    m_tree->Branch("posy_mult", &posy_mult);
    m_tree->Branch("posz_mult", &posz_mult);
    m_tree->Branch("e_mult", &e_mult);
    m_tree->Branch("sec_mult", &sec_mult);
    m_tree->Branch("nsec_mult", &nsec_mult);
    m_tree->Branch("mult_type", &mult_type);

    m_tree->Branch("totq_mult", &totq_mult);
    m_tree->Branch("meanq_mult", &meanq_mult);
    m_tree->Branch("stdq_mult", &stdq_mult);
    m_tree->Branch("minq_mult", &minq_mult);
    m_tree->Branch("maxq_mult", &maxq_mult);
    m_tree->Branch("nhit_mult", &nhit_mult);
    m_tree->Branch("meant_mult", &meant_mult);
    m_tree->Branch("stdt_mult", &stdt_mult);

    m_tree->Branch("method_mu", &method_mu);
    m_tree->Branch("loc_mu", &loc_mu);
    m_tree->Branch("posx_mu", &posx_mu);
    m_tree->Branch("posy_mu", &posy_mu);
    m_tree->Branch("posz_mu", &posz_mu);
    m_tree->Branch("dirx_mu", &dirx_mu);
    m_tree->Branch("diry_mu", &diry_mu);
    m_tree->Branch("dirz_mu", &dirz_mu);
    m_tree->Branch("totq_mu", &totq_mu);
    m_tree->Branch("sec_mu", &sec_mu);
    m_tree->Branch("nsec_mu", &nsec_mu);
    m_tree->Branch("quality_mu", &quality_mu);

    return true;
}

void IBDAnalysis::process(const EventContext::View& events) {
    run_id = events.runid();
    std::vector<TimeRangeMuonVetoSelection> mu_cut;
    std::vector<TimeRangeMuonVetoSelection> mu_spa_neu_cut;
    for (const track& trk : events.tracks()) {
        mu_cut.emplace_back(trk, TimeStamp{0, 0}, TimeStamp{0, 5000000});
        mu_spa_neu_cut.emplace_back(trk, TimeStamp{0, 20000}, TimeStamp{0, 2000000});
    }

    FiducialVolumeSelection fiducial_vol_cut{17700.0 /* 16500.0 */};
    // ChimneySelection chimney_cut{15500.0, 3000.0};
    // ChargeRangeSelection prompt_charge_cut{1500.0, 20000.0};
    // ChargeRangeSelection delayed_charge_cut{4000.0, 6000.0};
    EnergyRangeSelection prompt_energy_cut {0.6, 20.0 /* 0.7, 12.0 */};
    EnergyRangeSelection delayed_energy_cut{0.6, 3.0 /* 2.0, 2.5 */};
    // EnergyRangeSelection spa_neu_energy_cut{1.5, 20.0};
    // EnergyRangeSelection multiplicity_energy_cut{2.0, 12.0};

    std::vector<VertexCorrelationSelection> spa_neu_cut;
    for (const vertex& neu : events.vertices()) {
        bool is_in_veto = false;
        for (const TimeRangeMuonVetoSelection& cut : mu_spa_neu_cut) {
            if (!cut.isIn(neu)) continue;
            is_in_veto = true;
            break;
        }
        if (!is_in_veto) continue;
        if (!prompt_energy_cut.isIn(neu)) continue;
        spa_neu_cut.emplace_back(neu, 40000.0, TimeStamp{-2000000000} /* TimeStamp{0, 0} */, TimeStamp{2000000000});
    }

    std::vector<ibd_info> ibds;

    for (const vertex& prompt : events.current()) {
        LogInfo << prompt << '\n';
        if (!fiducial_vol_cut.isIn(prompt)) {
            LogInfo << "Prompt not in fiducial volume\n";
            continue;
        }
        // if (chimney_cut.isIn(prompt)) {
        //     LogInfo << "Prompt is a chimney\n";
        //     continue;
        // }

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

            ibd_info cand(prompt, delayed);
            LogInfo << "IBD event detected!\n";

            for (const VertexCorrelationSelection& cut : spa_neu_cut) {
                if (!cut.isIn(cand.pair.prompt) || !cut.isIn(cand.pair.delayed)) continue;
                cand.neus.push_back(cut.c_vtx);
            }

            TimeRangeSelection multi_prompt_time{prompt.ts, TimeStamp{0, -1000000}, TimeStamp{0, 0}};
            TimeRangeSelection multi_between_time{prompt.ts, TimeStamp{0, 0}, delayed.ts - prompt.ts};
            TimeRangeSelection multi_delayed_time{delayed.ts, TimeStamp{0, 0}, TimeStamp{0, 1000000}};
            for (const vertex& mult : events.vertices()) {
                if (mult.ts == cand.pair.prompt.ts) continue;
                if (mult.ts == cand.pair.delayed.ts) continue;

                // if (!fiducial_vol_cut.isIn(mult)) continue;
                // if (chimney_cut.isIn(mult)) continue;

                if (!prompt_energy_cut.isIn(mult)) continue;

                is_vetoed = false;
                for (const TimeRangeMuonVetoSelection& cut : mu_cut) {
                    if (!cut.isIn(mult)) continue;
                    is_vetoed = true;
                    break;
                }
                if (is_vetoed) continue;

                if (multi_prompt_time.isIn(mult)) {
                    cand.mults.push_back({mult, 0});
                }
                else if (multi_between_time.isIn(mult)) {
                    cand.mults.push_back({mult, 1});
                }
                else if (multi_delayed_time.isIn(mult)) {
                    cand.mults.push_back({mult, 2});
                }
                else {

                }
            }

            ibds.push_back(std::move(cand));
        }
    }

    method_mu.clear();
    loc_mu.clear();
    posx_mu.clear();
    posy_mu.clear();
    posz_mu.clear();
    dirx_mu.clear();
    diry_mu.clear();
    dirz_mu.clear();
    totq_mu.clear();
    sec_mu.clear();
    nsec_mu.clear();
    quality_mu.clear();

    for (const track& trk : events.tracks()) {
        vec3 dir = unit(trk.fpos - trk.ipos);
        method_mu.push_back(trk.method);
        loc_mu.push_back(trk.det);
        posx_mu.push_back(trk.ipos.x);
        posy_mu.push_back(trk.ipos.y);
        posz_mu.push_back(trk.ipos.z);
        dirx_mu.push_back(dir.x);
        diry_mu.push_back(dir.y);
        dirz_mu.push_back(dir.z);
        totq_mu.push_back(trk.totpe);
        sec_mu.push_back(trk.ts.GetSec());
        nsec_mu.push_back(trk.ts.GetNanoSec());
        quality_mu.push_back(trk.quality);
    }

    for (const ibd_info& ibd : ibds) {
        posx_p = ibd.pair.prompt.pos.x;
        posy_p = ibd.pair.prompt.pos.y;
        posz_p = ibd.pair.prompt.pos.z;
        e_p = ibd.pair.prompt.energy;
        sec_p = ibd.pair.prompt.ts.GetSec();
        nsec_p = ibd.pair.prompt.ts.GetNanoSec();

        totq_p = ibd.pair.prompt.calib.totq;
        meanq_p = ibd.pair.prompt.calib.meanq;
        stdq_p = ibd.pair.prompt.calib.stdq;
        minq_p = ibd.pair.prompt.calib.minq;
        maxq_p = ibd.pair.prompt.calib.maxq;
        nhit_p = ibd.pair.prompt.calib.nhit;
        meant_p = ibd.pair.prompt.calib.meant;
        stdt_p = ibd.pair.prompt.calib.stdt;

        posx_d = ibd.pair.delayed.pos.x;
        posy_d = ibd.pair.delayed.pos.y;
        posz_d = ibd.pair.delayed.pos.z;
        e_d = ibd.pair.delayed.energy;
        sec_d = ibd.pair.delayed.ts.GetSec();
        nsec_d = ibd.pair.delayed.ts.GetNanoSec();

        totq_d = ibd.pair.delayed.calib.totq;
        meanq_d = ibd.pair.delayed.calib.meanq;
        stdq_d = ibd.pair.delayed.calib.stdq;
        minq_d = ibd.pair.delayed.calib.minq;
        maxq_d = ibd.pair.delayed.calib.maxq;
        nhit_d = ibd.pair.delayed.calib.nhit;
        meant_d = ibd.pair.delayed.calib.meant;
        stdt_d = ibd.pair.delayed.calib.stdt;

        posx_n.clear();
        posy_n.clear();
        posz_n.clear();
        e_n.clear();
        totq_n.clear();
        sec_n.clear();
        nsec_n.clear();

        for (const vertex& neu : ibd.neus) {
            posx_n.push_back(neu.pos.x);
            posy_n.push_back(neu.pos.y);
            posz_n.push_back(neu.pos.z);
            e_n.push_back(neu.energy);
            sec_n.push_back(neu.ts.GetSec());
            nsec_n.push_back(neu.ts.GetNanoSec());

            totq_n.push_back(neu.calib.totq);
            meanq_n.push_back(neu.calib.meanq);
            stdq_n.push_back(neu.calib.stdq);
            minq_n.push_back(neu.calib.minq);
            maxq_n.push_back(neu.calib.maxq);
            nhit_n.push_back(neu.calib.nhit);
            meant_n.push_back(neu.calib.meant);
            stdt_n.push_back(neu.calib.stdt);
        }

        posx_mult.clear();
        posy_mult.clear();
        posz_mult.clear();
        e_mult.clear();
        totq_mult.clear();
        sec_mult.clear();
        nsec_mult.clear();
        mult_type.clear();

        for (const mult_info& mult : ibd.mults) {
            posx_mult.push_back(mult.vtx.pos.x);
            posy_mult.push_back(mult.vtx.pos.y);
            posz_mult.push_back(mult.vtx.pos.z);
            e_mult.push_back(mult.vtx.energy);
            sec_mult.push_back(mult.vtx.ts.GetSec());
            nsec_mult.push_back(mult.vtx.ts.GetNanoSec());
            mult_type.push_back(mult.type);

            totq_mult.push_back(mult.vtx.calib.totq);
            meanq_mult.push_back(mult.vtx.calib.meanq);
            stdq_mult.push_back(mult.vtx.calib.stdq);
            minq_mult.push_back(mult.vtx.calib.minq);
            maxq_mult.push_back(mult.vtx.calib.maxq);
            nhit_mult.push_back(mult.vtx.calib.nhit);
            meant_mult.push_back(mult.vtx.calib.meant);
            stdt_mult.push_back(mult.vtx.calib.stdt);
        }

        m_tree->Fill();
    }
}