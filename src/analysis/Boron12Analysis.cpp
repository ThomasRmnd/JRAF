#include "analysis/Boron12Analysis.hpp"

#include "SniperKernel/SniperLog.h"

#include "selection/Energy.hpp"
#include "selection/Muon.hpp"

Boron12Analysis::Boron12Analysis(const std::string& name, const std::string& method) : 
    Analysis{name, method} 
{}

bool Boron12Analysis::initialize() {
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

void Boron12Analysis::process(const EventContext::View& events) {
    run_id = events.runid();
    std::vector<TimeRangeMuonVetoSelection> mu_spa_neu_cut;
    for (const track& trk : events.tracks()) {
        mu_spa_neu_cut.emplace_back(trk, TimeStamp{0, 20000}, TimeStamp{0, 2000000});
    }

    FiducialVolumeSelection fiducial_vol_cut{18000.0};
    EnergyRangeSelection prompt_energy_cut {0.6, 20.0};

    for (const vertex& b12 : events.current()) {
        if (!fiducial_vol_cut.isIn(b12)) continue;

        if (!prompt_energy_cut.isIn(b12)) continue;
        bool is_in_veto = false;
        for (const TimeRangeMuonVetoSelection& cut : mu_spa_neu_cut) {
            if (!cut.isIn(b12)) continue;
            is_in_veto = true;
            break;
        }
        if (!is_in_veto) continue;

        posx = b12.pos.x;
        posy = b12.pos.y;
        posz = b12.pos.z;
        e = b12.energy;
        sec = b12.ts.GetSec();
        nsec = b12.ts.GetNanoSec();

        totq = b12.calib.totq;
        meanq = b12.calib.meanq;
        stdq = b12.calib.stdq;
        minq = b12.calib.minq;
        maxq = b12.calib.maxq;
        meant = b12.calib.meant;
        stdt = b12.calib.stdt;
        npmt = b12.calib.npmt;
        nhit = b12.calib.nhit;
        meanhit = b12.calib.meanhit;
        stdhit = b12.calib.stdhit;

        m_tree->Fill();
    }
}