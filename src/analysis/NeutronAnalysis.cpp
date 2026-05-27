#include "analysis/NeutronAnalysis.hpp"

#include "SniperKernel/SniperLog.h"

#include "selection/Energy.hpp"
#include "selection/Muon.hpp"

NeutronAnalysis::NeutronAnalysis(const std::string& name, const std::string& method) : 
    Analysis{name, method} 
{}

bool NeutronAnalysis::initialize() {
    if (!Analysis::initialize()) return false;

    m_tree->Branch("run_id", &run_id);

    m_tree->Branch("posx", &posx);
    m_tree->Branch("posy", &posy);
    m_tree->Branch("posz", &posz);
    m_tree->Branch("e", &e);
    m_tree->Branch("sec", &sec);
    m_tree->Branch("nsec", &nsec);

    m_tree->Branch("totq", &totq);
    m_tree->Branch("meanq", &meanq);
    m_tree->Branch("stdq", &stdq);
    m_tree->Branch("minq", &minq);
    m_tree->Branch("maxq", &maxq);
    m_tree->Branch("meant", &meant);
    m_tree->Branch("stdt", &stdt);
    m_tree->Branch("npmt", &npmt);
    m_tree->Branch("nhit", &nhit);
    m_tree->Branch("meanhit", &meanhit);
    m_tree->Branch("stdhit", &stdhit);

    return true;
}

void NeutronAnalysis::process(const EventContext::View& events) {
    run_id = events.runid();
    std::vector<TimeRangeMuonVetoSelection> mu_spa_neu_cut;
    for (const track& trk : events.tracks()) {
        mu_spa_neu_cut.emplace_back(trk, TimeStamp{0, 20000}, TimeStamp{0, 2000000});
    }

    FiducialVolumeSelection fiducial_vol_cut{18000.0};
    EnergyRangeSelection prompt_energy_cut {0.6, 20.0};

    for (const vertex& neu : events.current()) {
        if (!fiducial_vol_cut.isIn(neu)) continue;

        if (!prompt_energy_cut.isIn(neu)) continue;
        bool is_in_veto = false;
        for (const TimeRangeMuonVetoSelection& cut : mu_spa_neu_cut) {
            if (!cut.isIn(neu)) continue;
            is_in_veto = true;
            break;
        }
        if (!is_in_veto) continue;

        posx = neu.pos.x;
        posy = neu.pos.y;
        posz = neu.pos.z;
        e = neu.energy;
        sec = neu.ts.GetSec();
        nsec = neu.ts.GetNanoSec();

        totq = neu.calib.totq;
        meanq = neu.calib.meanq;
        stdq = neu.calib.stdq;
        minq = neu.calib.minq;
        maxq = neu.calib.maxq;
        meant = neu.calib.meant;
        stdt = neu.calib.stdt;
        npmt = neu.calib.npmt;
        nhit = neu.calib.nhit;
        meanhit = neu.calib.meanhit;
        stdhit = neu.calib.stdhit;

        m_tree->Fill();
    }
}