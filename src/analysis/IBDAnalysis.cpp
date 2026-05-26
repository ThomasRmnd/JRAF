#include "analysis/IBDAnalysis.hpp"

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
    m_tree->Branch("meant_n", &meant_n);
    m_tree->Branch("stdt_n", &stdt_n);
    m_tree->Branch("npmt_n", &npmt_n);
    m_tree->Branch("nhit_n", &nhit_n);
    m_tree->Branch("meanhit_n", &meanhit_n);
    m_tree->Branch("stdhit_n", &stdhit_n);

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
    m_tree->Branch("meant_mult", &meant_mult);
    m_tree->Branch("stdt_mult", &stdt_mult);
    m_tree->Branch("npmt_mult", &npmt_mult);
    m_tree->Branch("nhit_mult", &nhit_mult);
    m_tree->Branch("meanhit_mult", &meanhit_mult);
    m_tree->Branch("stdhit_mult", &stdhit_mult);

    m_tree->Branch("method_mu", &method_mu);
    m_tree->Branch("loc_mu", &loc_mu);
    m_tree->Branch("iposx_mu", &iposx_mu);
    m_tree->Branch("iposy_mu", &iposy_mu);
    m_tree->Branch("iposz_mu", &iposz_mu);
    m_tree->Branch("fposx_mu", &fposx_mu);
    m_tree->Branch("fposy_mu", &fposy_mu);
    m_tree->Branch("fposz_mu", &fposz_mu);
    m_tree->Branch("totq_cd_mu", &totq_cd_mu);
    m_tree->Branch("totq_wp_mu", &totq_wp_mu);
    m_tree->Branch("sec_mu", &sec_mu);
    m_tree->Branch("nsec_mu", &nsec_mu);
    m_tree->Branch("quality_mu", &quality_mu);

    m_tree_cutflow = new TTree((m_name + "__CutFlow").c_str(), (m_name + "__CutFlow").c_str());
    if (!m_tree_cutflow) {
        std::cout << "[ERROR] Could not create the TTree for CutFlow of Analysis " << m_name << '\n';
        return false;
    }

    m_tree_cutflow->Branch("run_id", &run_id);
    m_tree_cutflow->Branch("cf_prompt_total", &cf_prompt_total);
    m_tree_cutflow->Branch("cf_prompt_fv", &cf_prompt_fv);
    m_tree_cutflow->Branch("cf_prompt_energy", &cf_prompt_energy);
    m_tree_cutflow->Branch("cf_prompt_muon", &cf_prompt_muon);
    m_tree_cutflow->Branch("cf_pair_total", &cf_pair_total);
    m_tree_cutflow->Branch("cf_pair_delayed_fv", &cf_pair_delayed_fv);
    m_tree_cutflow->Branch("cf_pair_delayed_energy", &cf_pair_delayed_energy);
    m_tree_cutflow->Branch("cf_pair_corr", &cf_pair_corr);
    m_tree_cutflow->Branch("cf_pair_delayed_muon", &cf_pair_delayed_muon);
    m_tree_cutflow->Branch("cf_ibd_final", &cf_ibd_final);

    return true;
}

void IBDAnalysis::process(const EventContext::View& events) {
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

    std::vector<ibd_info> ibds;

    for (const vertex& prompt : events.current()) {
        // LogInfo << prompt << '\n';
        ++cf_prompt_total;
        
        if (!fiducial_vol_cut.isIn(prompt)) continue;
        ++cf_prompt_fv;

        if (!prompt_energy_cut.isIn(prompt)) continue;
        ++cf_prompt_energy;

        bool is_vetoed = false;
        for (const TimeRangeMuonVetoSelection& cut : mu_cut) {
            if (!cut.isIn(prompt)) continue;
            is_vetoed = true;
            break;
        }
        if (is_vetoed) continue;

        ++cf_prompt_muon;

        VertexCorrelationSelection correlation_cut{prompt, 1500.0, TimeStamp{0, 1000}, TimeStamp{0, 1500000}};

        for (const vertex& delayed : events.after()) {
            // LogInfo << delayed << '\n';
            ++cf_pair_total;

            if (!correlation_cut.isIn(delayed)) continue;
            ++cf_pair_corr;

            if (!fiducial_vol_cut.isIn(delayed)) continue;
            ++cf_pair_delayed_fv;

            if (!delayed_energy_cut_hydrogen.isIn(delayed) && !delayed_energy_cut_carbon.isIn(delayed)) continue;

            ++cf_pair_delayed_energy;

            is_vetoed = false;
            for (const TimeRangeMuonVetoSelection& cut : mu_cut) {
                if (!cut.isIn(delayed)) continue;
                is_vetoed = true;
                break;
            }
            if (is_vetoed) continue;
            ++cf_pair_delayed_muon;

            ++cf_ibd_final;

            ibd_info cand(prompt, delayed);
            LogInfo << "IBD event detected!\n";
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

    for (ibd_info& cand : ibds) {
        for (const VertexCorrelationSelection& cut : spa_neu_cut) {
            if (!cut.isIn(cand.pair.prompt) && !cut.isIn(cand.pair.delayed)) continue;
            cand.neus.push_back(cut.c_vtx);
        }

        TimeRangeSelection multi_prompt_time{cand.pair.prompt.ts, TimeStamp{0, -1000000}, TimeStamp{0, 0}};
        TimeRangeSelection multi_between_time{cand.pair.prompt.ts, TimeStamp{0, 0}, cand.pair.delayed.ts - cand.pair.prompt.ts};
        TimeRangeSelection multi_delayed_time{cand.pair.delayed.ts, TimeStamp{0, 0}, TimeStamp{0, 1000000}};
        for (const vertex& mult : events.vertices()) {
            if (mult.ts == cand.pair.prompt.ts) continue;
            if (mult.ts == cand.pair.delayed.ts) continue;

            // if (!fiducial_vol_cut.isIn(mult)) continue;
            // if (chimney_cut.isIn(mult)) continue;

            if (!prompt_energy_cut.isIn(mult)) continue;

            bool is_vetoed = false;
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
    }

    method_mu.clear();
    loc_mu.clear();
    iposx_mu.clear();
    iposy_mu.clear();
    iposz_mu.clear();
    fposx_mu.clear();
    fposy_mu.clear();
    fposz_mu.clear();
    totq_cd_mu.clear();
    totq_wp_mu.clear();
    sec_mu.clear();
    nsec_mu.clear();
    quality_mu.clear();

    for (const track& trk : events.tracks()) {
        method_mu.push_back(trk.method);
        loc_mu.push_back(trk.det);
        iposx_mu.push_back(trk.ipos.x);
        iposy_mu.push_back(trk.ipos.y);
        iposz_mu.push_back(trk.ipos.z);
        fposx_mu.push_back(trk.fpos.x);
        fposy_mu.push_back(trk.fpos.y);
        fposz_mu.push_back(trk.fpos.z);
        totq_cd_mu.push_back(trk.totq_cd);
        totq_wp_mu.push_back(trk.totq_wp);
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
        meant_p = ibd.pair.prompt.calib.meant;
        stdt_p = ibd.pair.prompt.calib.stdt;
        npmt_p = ibd.pair.prompt.calib.npmt;
        nhit_p = ibd.pair.prompt.calib.nhit;
        meanhit_p = ibd.pair.prompt.calib.meanhit;
        stdhit_p = ibd.pair.prompt.calib.stdhit;

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
        meant_d = ibd.pair.delayed.calib.meant;
        stdt_d = ibd.pair.delayed.calib.stdt;
        npmt_d = ibd.pair.delayed.calib.npmt;
        nhit_d = ibd.pair.delayed.calib.nhit;
        meanhit_d = ibd.pair.delayed.calib.meanhit;
        stdhit_d = ibd.pair.delayed.calib.stdhit;

        posx_n.clear();
        posy_n.clear();
        posz_n.clear();
        e_n.clear();
        sec_n.clear();
        nsec_n.clear();

        totq_n.clear();
        meanq_n.clear();
        stdq_n.clear();
        minq_n.clear();
        maxq_n.clear();
        meant_n.clear();
        stdt_n.clear();
        npmt_n.clear();
        nhit_n.clear();
        meanhit_n.clear();
        stdhit_n.clear();

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
            meant_n.push_back(neu.calib.meant);
            stdt_n.push_back(neu.calib.stdt);
            npmt_n.push_back(neu.calib.npmt);
            nhit_n.push_back(neu.calib.nhit);
            meanhit_n.push_back(neu.calib.meanhit);
            stdhit_n.push_back(neu.calib.stdhit);
        }

        posx_mult.clear();
        posy_mult.clear();
        posz_mult.clear();
        e_mult.clear();
        sec_mult.clear();
        nsec_mult.clear();
        mult_type.clear();

        totq_mult.clear();
        meanq_mult.clear();
        stdq_mult.clear();
        minq_mult.clear();
        maxq_mult.clear();
        meant_mult.clear();
        stdt_mult.clear();
        npmt_mult.clear();
        nhit_mult.clear();
        meanhit_mult.clear();
        stdhit_mult.clear();

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
            meant_mult.push_back(mult.vtx.calib.meant);
            stdt_mult.push_back(mult.vtx.calib.stdt);
            npmt_mult.push_back(mult.vtx.calib.npmt);
            nhit_mult.push_back(mult.vtx.calib.nhit);
            meanhit_mult.push_back(mult.vtx.calib.meanhit);
            stdhit_mult.push_back(mult.vtx.calib.stdhit);
        }

        m_tree->Fill();
    }
}

bool IBDAnalysis::write() {
    if (!Analysis::write()) return false;

    if (m_tree_cutflow) {
        m_tree_cutflow->Fill();
        m_tree_cutflow->Write();
    }

    return true;
}