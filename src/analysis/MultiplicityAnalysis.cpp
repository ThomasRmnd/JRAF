#include "analysis/MultiplicityAnalysis.hpp"

#include "SniperKernel/SniperLog.h"

#include "selection/Energy.hpp"
#include "selection/Muon.hpp"

MultiplicityAnalysis::MultiplicityAnalysis(const std::string& name, const std::string& method) : 
    Analysis{name, method}
{}

bool MultiplicityAnalysis::initialize() {
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

void MultiplicityAnalysis::process(const EventContext::View& events) {
    run_id = events.runid();
    std::vector<TimeRangeMuonVetoSelection> mu_cut;
    for (const track& trk : events.tracks()) {
        mu_cut.emplace_back(trk, TimeStamp{0, 0}, TimeStamp{0, 5000000});
    }

    FiducialVolumeSelection fiducial_vol_cut{18000.0};
    EnergyRangeSelection multiplicity_energy_cut {2.0, 12.0};

    for (const vertex& multi : events.current()) {
        if (!fiducial_vol_cut.isIn(multi)) continue;

        if (!multiplicity_energy_cut.isIn(multi)) continue;
        bool is_in_veto = false;
        for (const TimeRangeMuonVetoSelection& cut : mu_cut) {
            if (!cut.isIn(multi)) continue;
            is_in_veto = true;
            break;
        }
        if (is_in_veto) continue;

        posx = multi.pos.x;
        posy = multi.pos.y;
        posz = multi.pos.z;
        e = multi.energy;
        sec = multi.ts.GetSec();
        nsec = multi.ts.GetNanoSec();
        totq = multi.calib.totq;
        meanq = multi.calib.meanq;
        stdq = multi.calib.stdq;
        minq = multi.calib.minq;
        maxq = multi.calib.maxq;
        meant = multi.calib.meant;
        stdt = multi.calib.stdt;
        npmt = multi.calib.npmt;
        nhit = multi.calib.nhit;
        meanhit = multi.calib.meanhit;
        stdhit = multi.calib.stdhit;

        m_tree->Fill();
    }
}