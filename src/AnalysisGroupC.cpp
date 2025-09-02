#include "AnalysisGroupC.hpp"

#include <cmath>
#include <numeric>

#include "SniperKernel/AlgFactory.h"
#include "SniperKernel/SniperLog.h"

#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdTrackRecHeader.h"
#include "Event/CdTriggerHeader.h"
#include "Event/WpCalibHeader.h"
#include "Event/WpTriggerHeader.h"
#include "EvtNavigator/EvtNavHelper.h"

#include "analysis/FirstCrossCheckAnalysis.hpp"
#include "loader/BasicLoader.hpp"
#include "loader/JointLoader.hpp"
#include "loader/cd/CdLRangeFiller.hpp"
#include "loader/cd/CdRangeFiller.hpp"
#include "loader/cd/CdSRangeFiller.hpp"
#include "loader/tt/TtRangeFiller.hpp"
#include "loader/wp/WpRangeFiller.hpp"

DECLARE_ALGORITHM(AnalysisGroupC);

AnalysisGroupC::AnalysisGroupC(const std::string& name) : 
    AlgBase{name},
    m_iEvt{0u}
{
    declProp("RecTool", m_recToolName);
    declProp("Pmt3inchTimeReso", m_sigmaPmt3inch = 1.0);
    declProp("Pmt20inchTimeReso", m_sigmaPmt20inch = 8.0);
    declProp("PmtTTTimeReso", m_sigmaPmtTt = 2.0); // The sigma is not true, a placeholder
    declProp("Use3inchPMT", m_flagUse3inch = true);
    declProp("Use20inchPMT", m_flagUse20inch = true);
    declProp("ChosenDetectors", m_chosenDetectors = 0b111);
    declProp("UseJointLoader", m_useJointLoader = false);
    declProp("LoaderTimeWindow", m_loaderTimeWindow = {-500.0, 500.0});
    declProp("ReconstructMuonMode", m_reconstruct_muon_mode = false);
    declProp("OutputFilename", m_ofilename = "output.root");
}

bool AnalysisGroupC::initialize() {
    m_params.set("Pmt20inchTimeReso", m_sigmaPmt20inch);
    m_params.set("Pmt3inchTimeReso", m_sigmaPmt3inch);
    m_params.set("PmtTTTimeReso", m_sigmaPmtTt);

    if (!initBufSvc()) return false;
    if (!initGeomSvc()) return false;
    if (!initLoader()) return false;

    if (m_reconstruct_muon_mode) {
        if (!initRecTool()) return false;
    }
    else {
        if (!initAnalyses()) return false;
    }

    LogInfo  << objName() << " initialized successfully\n"; 
    return true;
}

bool AnalysisGroupC::initBufSvc() {
    SniperDataPtr<JM::NavBuffer> navBuf(getParent(), "/Event");
    if (navBuf.invalid()) {
        LogError << "Cannot get the NavBuffer @ /Event\n";
        return false;
    }
    m_buf = navBuf.data();
    return true; 
}

bool AnalysisGroupC::initGeomSvc() {
    SniperPtr<IRecGeomSvc> rgSvc(getParent(), "RecGeomSvc");
    if (rgSvc.invalid()) {
        LogError << "Failed to get RecGeomSvc instance!\n";
        return false;
    }
    m_rgSvc = rgSvc.data();

    SniperPtr<ITTGeomSvc> ttgSvc(getParent(), "TTGeomSvc");
    if (ttgSvc.invalid()) {
        LogError << "Cannot get the TTGeomSvc\n";
        return false;
    }
    m_ttgSvc = ttgSvc.data();
    return true;
}

bool AnalysisGroupC::initLoader() {
    std::shared_ptr<RangeFiller<CdGeom>> cd_filler = nullptr;
    std::shared_ptr<RangeFiller<WpGeom>> wp_filler = nullptr;
    std::shared_ptr<RangeFiller<TtGeom>> tt_filler = nullptr;

    DetectorType chosen_detector_type = static_cast<DetectorType>(m_chosenDetectors);

    if ( (chosen_detector_type & DetectorType::CD) == DetectorType::CD ) {
        if (m_flagUse20inch && m_flagUse3inch) 
            cd_filler = std::make_shared<CdRangeFiller>("CdRangeFiller", m_sigmaPmt20inch, m_sigmaPmt3inch);
        else if (m_flagUse20inch)
            cd_filler = std::make_shared<CdLRangeFiller>("CdLRangeFiller", m_sigmaPmt20inch);
        else if (m_flagUse3inch) 
            cd_filler = std::make_shared<CdSRangeFiller>("CdSRangeFiller", m_sigmaPmt3inch);
        // else cd_filler = nullptr;
    }
    if ( (chosen_detector_type & DetectorType::WP) == DetectorType::WP ) {
        wp_filler = std::make_shared<WpRangeFiller>("WpRangeFiller", m_sigmaPmt20inch);
    }
    if ( (chosen_detector_type & DetectorType::TT) == DetectorType::TT ) {
        tt_filler = std::make_shared<TtRangeFiller>("TtRangeFiller", m_sigmaPmtTt, m_ttgSvc);
    }

    if (m_useJointLoader)
        m_loader = std::make_shared<JointLoader>("JointLoader", &m_pmtTable, m_loaderTimeWindow, cd_filler, wp_filler, tt_filler, m_rgSvc);
    else
        m_loader = std::make_shared<BasicLoader>("BasicLoader", &m_pmtTable, cd_filler, wp_filler, tt_filler, m_rgSvc);

    if (!m_loader) {
        LogError << "Failed to create loader\n";
        return false;
    }
    
    return m_loader->initialize();
}

bool AnalysisGroupC::initRecTool() {
    m_recTool = tool<IRecMuonTool>(m_recToolName);
    if (!m_recTool) {
        LogError << "Failed to retrieve reconstruction tool named " << m_recToolName << '\n';
        return false;
    }
    if (!m_recTool->configure(&m_params, &m_pmtTable)) return false;
    // if (!m_recTool->initialize()) return false;
    return true;
}

bool AnalysisGroupC::initAnalyses() {
    m_file = TFile::Open(m_ofilename.c_str(), "RECREATE");
    if (!m_file) {
        LogError << "Failed to create output file: " << m_ofilename << '\n';
        return false;
    }
    m_analyses.push_back(std::make_shared<FirstCrossCheckAnalysis>("FirstCrossCheckAnalysis"));
    for (std::shared_ptr<Analysis>& analysis : m_analyses) {
        if (!analysis->initialize()) return false;
    }
    return true;
}

// Cuts for cross-check (as of August 29th)
// Muon Veto (CD+WP)	CD npe > 30k, WP npe > 400 with ±500 ns coincidence → 2 ms veto (50 μs dead time for CD, 4 μs dead time for WP)
// Muon Veto (WP-only)	WP npe > 400 → 2 ms veto (4 μs dead time)
// Muon Veto (CD-only)	CD npe > 30k, no WP trigger > 400 npe in ±500 ns and after >2 ms from the last tagged muon → 2 ms veto
// Job/files veto	Veto the first 2 ms for each job
// No header veto	2 ms veto after events without (OecHeader) & (CdLpmtCalibHeader or WpCalibHeader) & (CdTriggerHeader or WpTriggerHeader)
// Delayed NPE	Qdelayed∈[3700,6000]Qdelayed​∈[3700,6000]
// Prompt NPE	Qprompt∈[1500,20000]Qprompt​∈[1500,20000]
// Distance p–d	Δr<1.5 mΔr<1.5m
// Time p–d	1 μs<Δt<2 ms1μs<Δt<2ms
// Multiplicity	No events in prompt energy within 2 ms before the prompt or after the delayed event
// Fiducial Volume	r<17.2r<17.2 m ; (z<15.5z<15.5 m & ρ>3.0ρ>3.0 m)

bool AnalysisGroupC::execute() {
    LogInfo << "---------- Processing event by AnalysisGroupC: " << ++m_iEvt << " ----------\n";

    if (!m_loader->load(m_buf)) return false;

    std::size_t n_cd_used = 0ul, n_wp_used = 0ul, n_tt_used = 0ul;
    for (PmtProp& pmt : m_pmtTable) {
        if (!pmt.used) continue;
        if (pmt.loc == 1) ++n_cd_used;
        else if (pmt.loc == 2) ++n_wp_used;
        else if (pmt.loc == 3) ++n_tt_used;
    }

    LogInfo << "Number of PMTs in CD used: " << n_cd_used << '\n';
    LogInfo << "Number of PMTs in WP used: " << n_wp_used << '\n';
    LogInfo << "Number of PMTs in TT used: " << n_tt_used << '\n';

    double totq_cd = std::accumulate(m_pmtTable.begin(), m_pmtTable.end(), 0.0, 
        [](double sum, const PmtProp& pmt) { return sum + ( (pmt.used && (pmt.type & PmtType::PMT_20INCH) == pmt.type) ? pmt.q : 0.0 ); } 
    );
    double totq_wp = std::accumulate(m_pmtTable.begin(), m_pmtTable.end(), 0.0, 
        [](double sum, const PmtProp& pmt) { return sum + ( (pmt.used && (pmt.type & PmtType::PMT_WP) == pmt.type) ? pmt.q : 0.0 ); } 
    );

    JM::EvtNavigator* nav = m_buf->curEvt();
    TimeStamp ts{nav->TimeStamp().GetTimeSpec()};

    LogInfo << "TimeStamp: " << ts << '\n';
    LogInfo << "TotQ: CD = " << totq_cd << ", WP = " << totq_wp << '\n';

    if (m_reconstruct_muon_mode) {

        bool is_possibly_cd_muon = false;
        bool is_possibly_wp_muon = false;

        JM::CdLpmtCalibHeader* cd_lpmt_calib_hdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
        JM::CdTriggerHeader* cd_trig_hdr = JM::getHeaderObject<JM::CdTriggerHeader>(nav);
        JM::WpCalibHeader* wp_calib_hdr = JM::getHeaderObject<JM::WpCalibHeader>(nav);
        JM::WpTriggerHeader* wp_trig_hdr = JM::getHeaderObject<JM::WpTriggerHeader>(nav);
        if ((!cd_lpmt_calib_hdr || !cd_trig_hdr) && (!wp_calib_hdr || !wp_trig_hdr)) {
            is_possibly_cd_muon = true;
        }
        else if (
            totq_cd >= m_cd_muon_totq_thold && 
            totq_wp >= m_wp_muon_totq_thold && 
            ts - m_cd_last_muon > m_cd_afterpulse_thold &&
            ts - m_wp_last_muon > m_wp_afterpulse_thold
        ) {
            m_cd_last_muon = ts;
            m_wp_last_muon = ts;
            is_possibly_cd_muon = true;
            is_possibly_wp_muon = true;
        }
        else if (
            n_cd_used == 0ul && 
            totq_wp >= m_wp_muon_totq_thold &&
            ts - m_wp_last_muon > m_wp_afterpulse_thold
        ) {
            m_wp_last_muon = ts;
            is_possibly_wp_muon = true;
        }
        else if (
            totq_cd >= m_cd_muon_totq_thold && 
            totq_wp < m_wp_muon_totq_thold &&
            (ts - m_cd_last_muon > TimeStamp{0, 2000000} || ts - m_wp_last_muon > TimeStamp{0, 2000000})
        ) {
            m_cd_last_muon = ts;
            is_possibly_cd_muon = true;
        }

        /* RecTrks* trks = new RecTrks();
        if (!m_recTool->reconstruct(trks)) {
            LogWarn << "Failed to execute the reconstruction tool\n";
            return true;
        }

        if (!trks->size()) {
            LogInfo << "No track reconstructed\n";
            return true;
        } */

        if (is_possibly_cd_muon) {
            JM::CdTrackRecHeader* cd_hdr = new JM::CdTrackRecHeader();
            JM::CdTrackRecEvt* cd_evt = new JM::CdTrackRecEvt();

            JM::RecTrack* trk = new JM::RecTrack(
                CLHEP::HepLorentzVector(0.0, 0.0, 20000.0, 0.0),
                CLHEP::HepLorentzVector(0.0, 0.0, -20000.0, 0.0)
            );
            trk->setQuality(1.0f);
            trk->setPESum(1.0f);

            cd_evt->addTrack(trk);
            cd_hdr->setEvent(cd_evt);
            JM::addHeaderObject(nav, cd_hdr);
        }
        else if (is_possibly_wp_muon) {
            JM::CdTrackRecHeader* wp_hdr = new JM::CdTrackRecHeader();
            JM::CdTrackRecEvt* wp_evt = new JM::CdTrackRecEvt();

            JM::RecTrack* trk = new JM::RecTrack(
                CLHEP::HepLorentzVector(0.0, 0.0, 20000.0, 0.0),
                CLHEP::HepLorentzVector(0.0, 0.0, -20000.0, 0.0)
            );
            trk->setQuality(1.0f);
            trk->setPESum(1.0f);

            wp_evt->addTrack(trk);
            wp_hdr->setEvent(wp_evt);
            JM::addHeaderObject(nav, wp_hdr);
        }
    
    }

    else {
        for (std::shared_ptr<Analysis>& analysis : m_analyses) {
            analysis->process(m_buf);
        }
    }

    return true;
}

bool AnalysisGroupC::finalize() {
    if (m_loader && !m_loader->finalize()) return false;
    
    if (m_reconstruct_muon_mode) {
        if (m_recTool && !(dynamic_cast<ToolBase*>(m_recTool))->finalize()) return false;
    }
    else {
        if (!m_file) return false;
        m_file->cd();
        for (std::shared_ptr<Analysis>& ana : m_analyses) {
            if (!ana->write()) return false;
        }
        m_file->Close();
    }

    LogInfo  << objName() << " finalized successfully\n"; 
    return true;
}