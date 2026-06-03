#include "analysis/AccidentalAnalysis.hpp"

#include "SniperKernel/SniperLog.h"

#include "event/IBD.hpp"
#include "selection/Energy.hpp"
#include "selection/Muon.hpp"
#include "selection/Vertex.hpp"
#include "selection/Volume.hpp"

AccidentalAnalysis::AccidentalAnalysis(const std::string& name, const std::string& method) : 
    Analysis{name, method} 
{}

bool AccidentalAnalysis::initialize() {
    if (!Analysis::initialize()) return false;

    m_tree->Branch("run_id", &run_id);

    m_tree->Branch("posx_p", &posx_p);
    m_tree->Branch("posy_p", &posy_p);
    m_tree->Branch("posz_p", &posz_p);
    m_tree->Branch("e_p", &e_p);
    m_tree->Branch("sec_p", &sec_p);
    m_tree->Branch("nsec_p", &nsec_p);

    m_tree->Branch("totq_p", &totq_p);
    m_tree->Branch("meanq_p", &meanq_p);
    m_tree->Branch("stdq_p", &stdq_p);
    m_tree->Branch("minq_p", &minq_p);
    m_tree->Branch("maxq_p", &maxq_p);
    m_tree->Branch("meant_p", &meant_p);
    m_tree->Branch("stdt_p", &stdt_p);
    m_tree->Branch("npmt_p", &npmt_p);
    m_tree->Branch("nhit_p", &nhit_p);
    m_tree->Branch("meanhit_p", &meanhit_p);
    m_tree->Branch("stdhit_p", &stdhit_p);

    m_tree->Branch("posx_d", &posx_d);
    m_tree->Branch("posy_d", &posy_d);
    m_tree->Branch("posz_d", &posz_d);
    m_tree->Branch("e_d", &e_d);
    m_tree->Branch("sec_d", &sec_d);
    m_tree->Branch("nsec_d", &nsec_d);

    m_tree->Branch("totq_d", &totq_d);
    m_tree->Branch("meanq_d", &meanq_d);
    m_tree->Branch("stdq_d", &stdq_d);
    m_tree->Branch("minq_d", &minq_d);
    m_tree->Branch("maxq_d", &maxq_d);
    m_tree->Branch("meant_d", &meant_d);
    m_tree->Branch("stdt_d", &stdt_d);
    m_tree->Branch("npmt_d", &npmt_d);
    m_tree->Branch("nhit_d", &nhit_d);
    m_tree->Branch("meanhit_d", &meanhit_d);
    m_tree->Branch("stdhit_d", &stdhit_d);

    return true;
}

void AccidentalAnalysis::process(const EventContext::View& events) {
    run_id = events.runid();
    std::vector<TimeRangeMuonVetoSelection> mu_cut;
    std::vector<TimeRangeMuonVetoSelection> mu_spa_neu_cut;
    mu_cut.reserve(events.tracks().size());
    mu_spa_neu_cut.reserve(events.tracks().size());
    for (const track& trk : events.tracks()) {
        mu_cut.emplace_back(trk, TimeStamp{0, 0}, TimeStamp{0, 5000000});
        mu_spa_neu_cut.emplace_back(trk, TimeStamp{0, 20000}, TimeStamp{0, 2000000});
    }

    FiducialVolumeSelection fiducial_vol_cut{18000.0};
    EnergyRangeSelection prompt_energy_cut {0.6, 20.0};
    EnergyRangeSelection delayed_energy_cut_hydrogen{0.6, 3.0};
    EnergyRangeSelection delayed_energy_cut_carbon{4.0, 6.0};
    // EnergyRangeSelection spa_neu_energy_cut{1.5, 20.0};
    // EnergyRangeSelection multiplicity_energy_cut{2.0, 12.0};

    std::vector<ibd> ibds;

    for (const vertex& prompt : events.current()) {
        if (!fiducial_vol_cut.isIn(prompt)) continue;

        if (!prompt_energy_cut.isIn(prompt)) continue;

        bool is_vetoed = false;
        for (const TimeRangeMuonVetoSelection& cut : mu_cut) {
            if (!cut.isIn(prompt)) continue;
            is_vetoed = true;
            break;
        }
        if (is_vetoed) continue;

        VertexCorrelationSelection correlation_cut{prompt, 1500.0, TimeStamp{2, 0}, TimeStamp{4, 0}};

        for (const vertex& delayed : events.after()) {
            if (!correlation_cut.isIn(delayed)) continue;

            if (!fiducial_vol_cut.isIn(delayed)) continue;

            if (!delayed_energy_cut_hydrogen.isIn(delayed) && !delayed_energy_cut_carbon.isIn(delayed)) continue;

            is_vetoed = false;
            for (const TimeRangeMuonVetoSelection& cut : mu_cut) {
                if (!cut.isIn(delayed)) continue;
                is_vetoed = true;
                break;
            }
            if (is_vetoed) continue;

            ibd cand{prompt, delayed};
            LogInfo << "Accidental event detected!\n";
            ibds.push_back(std::move(cand));
        }
    }

    if (ibds.empty()) return;

    std::vector<VertexCorrelationSelection> spa_neu_cut;
    for (const vertex& neu : events.vertices()) {
        if (!prompt_energy_cut.isIn(neu)) continue;
        bool is_in_veto = false;
        for (const TimeRangeMuonVetoSelection& cut : mu_spa_neu_cut) {
            if (!cut.isIn(neu)) continue;
            is_in_veto = true;
            break;
        }
        if (!is_in_veto) continue;
        spa_neu_cut.emplace_back(neu, 40000.0, TimeStamp{0, -2000000000} /* TimeStamp{0, 0} */, TimeStamp{0, 2000000000});
    }

    for (const ibd& ibd : ibds) {
        posx_p = ibd.prompt.pos.x;
        posy_p = ibd.prompt.pos.y;
        posz_p = ibd.prompt.pos.z;
        e_p = ibd.prompt.energy;
        sec_p = ibd.prompt.ts.GetSec();
        nsec_p = ibd.prompt.ts.GetNanoSec();

        totq_p = ibd.prompt.calib.totq;
        meanq_p = ibd.prompt.calib.meanq;
        stdq_p = ibd.prompt.calib.stdq;
        minq_p = ibd.prompt.calib.minq;
        maxq_p = ibd.prompt.calib.maxq;
        meant_p = ibd.prompt.calib.meant;
        stdt_p = ibd.prompt.calib.stdt;
        npmt_p = ibd.prompt.calib.npmt;
        nhit_p = ibd.prompt.calib.nhit;
        meanhit_p = ibd.prompt.calib.meanhit;
        stdhit_p = ibd.prompt.calib.stdhit;

        posx_d = ibd.delayed.pos.x;
        posy_d = ibd.delayed.pos.y;
        posz_d = ibd.delayed.pos.z;
        e_d = ibd.delayed.energy;
        sec_d = ibd.delayed.ts.GetSec();
        nsec_d = ibd.delayed.ts.GetNanoSec();

        totq_d = ibd.delayed.calib.totq;
        meanq_d = ibd.delayed.calib.meanq;
        stdq_d = ibd.delayed.calib.stdq;
        minq_d = ibd.delayed.calib.minq;
        maxq_d = ibd.delayed.calib.maxq;
        meant_d = ibd.delayed.calib.meant;
        stdt_d = ibd.delayed.calib.stdt;
        npmt_d = ibd.delayed.calib.npmt;
        nhit_d = ibd.delayed.calib.nhit;
        meanhit_d = ibd.delayed.calib.meanhit;
        stdhit_d = ibd.delayed.calib.stdhit;

        m_tree->Fill();
    }
}

bool AccidentalAnalysis::write() {
    if (!Analysis::write()) return false;
    return true;
}