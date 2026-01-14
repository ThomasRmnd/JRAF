#ifndef ANALYSISGROUPC_ANALYSISGROUPC_HPP_
#define ANALYSISGROUPC_ANALYSISGROUPC_HPP_

#include "SniperKernel/AlgBase.h"

#include <deque>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TTree.h>

#include "Context/TimeStamp.h"
#include "Event/CdTrackRecHeader.h"
#include "Event/CdVertexRecHeader.h"
#include "Event/OecHeader.h"
#include "Event/TtRecHeader.h"
#include "Event/WpRecHeader.h"
#include "EvtNavigator/NavBuffer.h"
#include "Geometry/IPMTParamSvc.h"
#include "Geometry/IRecGeomSvc.hh"
#include "Geometry/ITTGeomSvc.hh"
#include "RootIOSvc/RootInputSvc.h"
#include "RecTools/IRecMuonTool.h"
#include "RecTools/PmtProp.h"
#include "UtilsThomas/loader/Loader.hpp"

#include "analysis/Analysis.hpp"
#include "event/Event.hpp"
#include "veto/Veto.hpp"

struct DAQTimeSaver {

    TTree* tree = nullptr;
    int run_id = 0;
    time_t sec = 0l;
    int nsec = 0;

    bool is_initialized = false;
    TimeStamp last_ts{0, 0};

    bool init() {
        tree = new TTree("DAQ", "DAQ");
        if (!tree) {
            LogError << "Cannot create DAQ TTree\n";
            return false;
        }
        tree->Branch("run_id", &run_id);
        tree->Branch("sec", &sec);
        tree->Branch("nsec", &nsec);
        return true;
    }

    bool add(const TimeStamp& ts, int run) {
        run_id = run;
        if (!is_initialized) {
            last_ts = ts;
            is_initialized = true;
            return true;
        }
        TimeStamp diff = ts - last_ts;
        TimeStamp daqtime{sec, nsec};
        daqtime.Add(diff);
        sec = daqtime.GetSec();
        nsec = daqtime.GetNanoSec();
        last_ts = ts;
        return true;
    };

    bool write() {
        if (!tree) return false;
        tree->Fill();
        tree->Write();
        return true;
    }

};

struct VetoTimeSaver {

    std::unordered_map<VetoType, TimeStamp> veto_map {
        {VetoType::BeginningOfJob, TimeStamp{1, 200000000}},
        {VetoType::MissingHeaders, TimeStamp{0, 5000000}},
        {VetoType::BigGaps, TimeStamp{1, 200000000}},
        {VetoType::Muon, TimeStamp{0, 5000000}}
    };

    TTree* tree = nullptr;

    int run_id = 0;
    time_t sec = 0l;
    int nsec = 0;

    unsigned char veto_type = 0;
    time_t veto_sec = 0l;
    int veto_nsec = 0;

    TimeStamp last_ts{0, 0};
    TimeStamp veto_duration{0, 0};

    bool init() {
        tree = new TTree("Veto", "Veto");
        if (!tree) {
            LogError << "Cannot create veto TTree\n";
            return false;
        }
        tree->Branch("run_id", &run_id);
        tree->Branch("sec", &sec);
        tree->Branch("nsec", &nsec);
        tree->Branch("veto_type", &veto_type);
        tree->Branch("veto_sec", &veto_sec);
        tree->Branch("veto_nsec", &veto_nsec);
        return true;
    }

    bool create(const TimeStamp& ts, VetoType type, int run) {
        std::unordered_map<VetoType, TimeStamp>::const_iterator it = veto_map.find(type);
        if (it == veto_map.end()) {
            return false;
        }
        run_id = run;
        sec = ts.GetSec();
        nsec = ts.GetNanoSec();
        last_ts = ts;
        veto_type = static_cast<unsigned char>(it->first);
        veto_duration = it->second;
        veto_sec = veto_duration.GetSec();
        veto_nsec = veto_duration.GetNanoSec();
        tree->Fill();
        return true;
    }

    bool create_no_veto(const TimeStamp& ts, VetoType type, int run) {
        std::unordered_map<VetoType, TimeStamp>::const_iterator it = veto_map.find(type);
        if (it == veto_map.end()) {
            return false;
        }
        run_id = run;
        sec = ts.GetSec();
        nsec = ts.GetNanoSec();
        veto_type = static_cast<unsigned char>(it->first);
        veto_sec = it->second.GetSec();
        veto_nsec = it->second.GetNanoSec();
        tree->Fill();
        return true;
    }

    bool inVeto(const TimeStamp& ts) {
        TimeStamp diff = ts - last_ts;
        if (diff < veto_duration) {
            return true;
        }
        return false;
    }

    bool write() {
        if (!tree) return false;
        tree->Write();
        return true;
    }

};

struct ContextFileTracker {

    std::string prevctx;
    std::string nextctx;

    std::string current;
    std::string next;
    bool change = false;

    bool isTarget(RootInputSvc* iptSvc) {
        if (current.empty()) {
            current = std::filesystem::path(iptSvc->getInputStream("EvtNavigator")->streamname()).filename().string();
            next = current;
        }
        if (change) {
            current = next;
            change = false;
        }
        std::string filename = std::filesystem::path(iptSvc->getInputStream("EvtNavigator")->streamname()).filename().string();
        if (filename != current) {
            next = filename;
            change = true;
        }
        LogInfo << "Current filename: " << current << '\n';
        LogInfo << "Previous filename: " << prevctx << '\n';
        LogInfo << "Next filename: " << nextctx << '\n';
        return (current != prevctx) && (current != nextctx);
    }

};

struct TtRecoFile {

    std::string filename;
    std::string treename = "TT";
    TChain* chain = nullptr;
    int cur_idx = 0;
    Int_t ntracks;
    TTimeStamp* ts = nullptr;
    Int_t npts[20];
    Double_t coeff0[20], coeff1[20], coeff2[20], coeff3[20], coeff4[20], coeff5[20];
    Double_t chi2[20];

    bool init() {
        chain = new TChain(treename.c_str());
        if (!chain) {
            LogError << "Cannot create chain with name " << treename << '\n';
            return false;
        }
        chain->Add(filename.c_str());

        chain->SetBranchAddress("NTracks", &ntracks);
        chain->SetBranchAddress("NPoints", &npts);
        chain->SetBranchAddress("start_TS", &ts);
        chain->SetBranchAddress("Coeff0", &coeff0);
        chain->SetBranchAddress("Coeff1", &coeff1);
        chain->SetBranchAddress("Coeff2", &coeff2);
        chain->SetBranchAddress("Coeff3", &coeff3);
        chain->SetBranchAddress("Coeff4", &coeff4);
        chain->SetBranchAddress("Coeff5", &coeff5);
        chain->SetBranchAddress("Chi2", &chi2);

        LogInfo << "TtRecoFile has " << chain->GetEntries() << " entries\n";

        return true;
    }

    bool find(const TimeStamp& ts_) {
        long nentries = chain->GetEntries();
        TimeStamp lower_bound = ts_ - TimeStamp{0, 1000};
        TimeStamp upper_bound = ts_;
        upper_bound.Add(TimeStamp{0, 1000});
        for (; cur_idx < nentries; ++cur_idx) {
            chain->GetEntry(cur_idx);
            
            TimeStamp cur_ts{ts->GetTimeSpec()};
            if (cur_ts < lower_bound) continue;
            else if (upper_bound < cur_ts) break;

            return true;
        }
        return false;
    }

};

struct TrackSaver {

    std::string filename;
    TFile* file = nullptr;
    TTree* tree = nullptr;

    int run_id = 0;
    time_t sec = 0l;
    int nsec = 0;

    std::vector<std::string> method;
    std::vector<unsigned char> det;
    std::vector<double> quality;

    std::vector<double> iposx;
    std::vector<double> iposy;
    std::vector<double> iposz;
    std::vector<double> fposx;
    std::vector<double> fposy;
    std::vector<double> fposz;


    bool init() {
        file = TFile::Open(filename.c_str(), "RECREATE");
        if (!file) {
            LogWarn << "Cannot open ROOT file " << filename << ". Skipping track saving\n";
            return true;
        }
        tree = new TTree("muons", "muons");
        tree->Branch("run_id", &run_id);
        tree->Branch("sec", &sec);
        tree->Branch("nsec", &nsec);
        tree->Branch("method", &method);
        tree->Branch("det", &det);
        tree->Branch("quality", &quality);
        tree->Branch("iposx", &iposx);
        tree->Branch("iposy", &iposy);
        tree->Branch("iposz", &iposz);
        tree->Branch("fposx", &fposx);
        tree->Branch("fposy", &fposy);
        tree->Branch("fposz", &fposz);
        return true;
    }

    void reset() {
        run_id = 0;
        sec = 0l;
        nsec = 0;
        method.clear();
        det.clear();
        quality.clear();
        iposx.clear();
        iposy.clear();
        iposz.clear();
        fposx.clear();
        fposy.clear();
        fposz.clear();
    }

    void add(RecTrks& trks_, const std::string& method_, int run_id_, const TimeStamp& ts_, const track::loc& force_loc = track::loc::none) {
        run_id = run_id_;
        sec = ts_.GetSec();
        nsec = ts_.GetNanoSec();
        for (int k = 0; k < trks_.size(); ++k) {
            method.push_back(method_);
            if (force_loc != track::loc::none) {
                det.push_back(static_cast<unsigned char>(force_loc));
            }
            else {
                det.push_back((trks_.isCdUsed(k) << 0) | (trks_.isWpUsed(k) << 1) | (trks_.isTtUsed(k) << 2));
            }
            quality.push_back(trks_.getQuality(k));
            iposx.push_back(trks_.getStart(k).X());
            iposy.push_back(trks_.getStart(k).Y());
            iposz.push_back(trks_.getStart(k).Z());
            fposx.push_back(trks_.getEnd(k).X());
            fposy.push_back(trks_.getEnd(k).Y());
            fposz.push_back(trks_.getEnd(k).Z());
        }
    }
    
    void add(JM::CdTrackRecHeader* cdt_hdr, const std::string& method_, int run_id_, const TimeStamp& ts_) {
        if (!cdt_hdr || !cdt_hdr->event()) return;
        const std::vector<JM::RecTrack*>& rec_tracks = cdt_hdr->event()->tracks();
        run_id = run_id_;
        sec = ts_.GetSec();
        nsec = ts_.GetNanoSec();
        for (JM::RecTrack* t : rec_tracks) {
            method.push_back(method_);
            det.push_back(0b001);
            quality.push_back(t->quality());
            iposx.push_back(t->start().x());
            iposy.push_back(t->start().y());
            iposz.push_back(t->start().z());
            fposx.push_back(t->end().x());
            fposy.push_back(t->end().y());
            fposz.push_back(t->end().z());
        }
    }

    void add(JM::WpRecHeader* wpt_hdr, const std::string& method_, int run_id_, const TimeStamp& ts_) {
        if (!wpt_hdr || !wpt_hdr->event()) return;
        const std::vector<JM::RecTrack*>& rec_tracks = wpt_hdr->event()->tracks();
        run_id = run_id_;
        sec = ts_.GetSec();
        nsec = ts_.GetNanoSec();
        for (JM::RecTrack* t : rec_tracks) {
            method.push_back(method_);
            det.push_back(0b010);
            quality.push_back(t->quality());
            iposx.push_back(t->start().x());
            iposy.push_back(t->start().y());
            iposz.push_back(t->start().z());
            fposx.push_back(t->end().x());
            fposy.push_back(t->end().y());
            fposz.push_back(t->end().z());
        }
    }

    void add(JM::TtRecHeader* ttt_hdr, const std::string& method_, int run_id_, const TimeStamp& ts_) {
        if (!ttt_hdr || !ttt_hdr->event()) return;
        JM::TtRecEvt* ttt_evt = ttt_hdr->event();
        for (int k = 0; k < ttt_evt->nTracks(); ++k) {
            vec3 ipos{ttt_evt->Coeff0()[k], ttt_evt->Coeff1()[k], ttt_evt->Coeff2()[k]};
            vec3 dir = unit(vec3{ttt_evt->Coeff3()[k], ttt_evt->Coeff4()[k], ttt_evt->Coeff5()[k]});
            vec3 fpos = ipos - 2.0 * dot(ipos, dir) * dir;
            method.push_back(method_);
            det.push_back(0b100);
            quality.push_back(ttt_evt->Chi2()[k]);
            iposx.push_back(ipos.x);
            iposy.push_back(ipos.y);
            iposz.push_back(ipos.z);
            fposx.push_back(fpos.x);
            fposy.push_back(fpos.y);
            fposz.push_back(fpos.z);
        }
    }

    void fill() {
        if (tree && !method.empty()) tree->Fill();
    }

    bool save() {
        if (!file || !tree) return true;
        file->cd();
        tree->Write();
        file->Close();
        return true;
    }

};

class NavBufferWrapper : public JM::NavBuffer {

public:

    NavBufferWrapper(JM::NavBuffer& buf) {
        m_dBuf.assign(buf.begin(), buf.end());
        m_iCur = 0;
    }

    ~NavBufferWrapper() override = default;

    void prev() {
        --m_iCur;
    }

    void next() {
        ++m_iCur;
    }

};

class AnalysisGroupC : public AlgBase {

public:

    AnalysisGroupC(const std::string& name);

    ~AnalysisGroupC() override = default;

    bool initialize() override;
    bool execute() override;
    bool finalize() override;

private:

    std::size_t m_iEvt = 0ul;
    TimeStamp m_tsEvt{0, 0};
    
    JM::NavBuffer* m_buf;
    RootInputSvc* m_iptSvc;

    // Reconstruction

    std::string m_recToolName;
	IRecMuonTool* m_recTool;
    std::string m_loaderName;
    std::string m_cdFillerName;
    std::string m_wpFillerName;
    std::string m_ttFillerName;
    Loader* m_loader;
    PmtTable m_pmtTable;
    Params m_params;

    // Muon selection variable

    double m_cd_muon_totq_thold = 30000.0;
    double m_wp_muon_totq_thold = 400.0;
    TimeStamp m_cd_afterpulse_thold{0, 50000};
    TimeStamp m_wp_afterpulse_thold{0, 4000};
    TimeStamp m_cd_last_muon{0, 0};
    TimeStamp m_wp_last_muon{0, 0};
    
    TtRecoFile m_ttRecoFile;
    ContextFileTracker m_contextTracker;
    BeginningOfJobVetoTracker m_begOfJobVetoTrkr;
    MissingHeaderVetoTracker m_missHdrVetoTrkr;
    BigGapsVetoTracker m_bigGapsVetoTrkr;
    MuonVetoTracker m_muvetoTrkr;

    std::string m_ofilename;
    TFile* m_file;
    DAQTimeSaver m_daqTimeSaver;
    VetoTimeSaver m_vetoTimeSaver;
    std::vector<std::string> m_methods;
    std::vector<std::shared_ptr<Analysis>> m_analyses;

    TrackSaver m_trkSaver;

	bool initBufSvc();
    bool initRecTool();
    bool initLoader();
    bool initAnalyses();

    void addTrack(RecTrks& rec_tracks, const std::string& method, const TimeStamp& ts, const track::loc& det, std::vector<track>& tracks);
    void addTrack(JM::CdTrackRecHeader* cdt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks);
    void addTrack(JM::WpRecHeader* wpt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks);
    void addTrack(JM::TtRecHeader* ttt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks);
    void addVertex(JM::OecHeader* oec_hdr, const std::string& method, const TimeStamp& ts, const calibration_context& calib, std::vector<vertex>& vertices);
    void addVertex(JM::CdVertexRecHeader* cdv_hdr, const std::string& method, const TimeStamp& ts, const calibration_context& calib, std::vector<vertex>& vertices);

};

#endif // ANALYSISGROUPC_ANALYSISGROUPC_HPP_