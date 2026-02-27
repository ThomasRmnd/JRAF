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
#include "Event/CalibPmtChannel.h"
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
#include "UtilsThomas/utils/DetectorType.hpp"

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

struct VetoWindow {
    TimeStamp start;
    TimeStamp end;
    VetoType type;
};

struct VetoTimeSaver {

    std::unordered_map<VetoType, TimeStamp> veto_map {
        {VetoType::BeginningOfJob, TimeStamp{1, 200000000}},
        {VetoType::MissingHeaders, TimeStamp{0, 5000000}},
        {VetoType::BigGaps, TimeStamp{1, 200000000}},
        {VetoType::Muon, TimeStamp{0, 5000000}}
    };

    std::deque<VetoWindow> active_vetoes;

    TTree* tree = nullptr;

    int run_id = 0;
    time_t sec = 0l;
    int nsec = 0;

    unsigned char veto_type = 0;
    time_t veto_sec = 0l;
    int veto_nsec = 0;

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

        TimeStamp duration = it->second;
        VetoWindow window{ts, ts + duration, type};
        active_vetoes.push_back(window);

        run_id = run;
        sec = ts.GetSec();
        nsec = ts.GetNanoSec();
        veto_type = static_cast<unsigned char>(it->first);
        veto_sec = duration.GetSec();
        veto_nsec = duration.GetNanoSec();
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
        cleanupExpired(ts);
        for (const VetoWindow& veto : active_vetoes) {
            if (veto.start <= ts && ts < veto.end) {
                return true;
            }
        }
        return false;
    }

    void cleanupExpired(const TimeStamp& ts) {
        while (!active_vetoes.empty() && ts >= active_vetoes.front().end) {
            active_vetoes.pop_front();
        }
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
    Long64_t cur_idx = 0;
    bool first_search = true;
    
    Int_t evtID, NTotPoints, NTracks;
    Int_t NPoints[20];
    TTimeStamp* start_TS = nullptr;
    Float_t PointX[20], PointY[20], PointZ[20];
    Double_t Coeff0[20], Coeff1[20], Coeff2[20], Coeff3[20], Coeff4[20], Coeff5[20];
    Double_t Chi2[20];

    bool init() {
        chain = new TChain(treename.c_str());
        if (!chain) {
            LogError << "Cannot create chain with name " << treename << '\n';
            return false;
        }
        chain->Add(filename.c_str());

        chain->SetBranchAddress("evtID", &evtID);
        chain->SetBranchAddress("NTotPoints", &NTotPoints);
        chain->SetBranchAddress("PointX", &PointX);
        chain->SetBranchAddress("PointY", &PointY);
        chain->SetBranchAddress("PointZ", &PointZ);
        chain->SetBranchAddress("NTracks", &NTracks);
        chain->SetBranchAddress("NPoints", NPoints);
        chain->SetBranchAddress("start_TS", &start_TS);
        chain->SetBranchAddress("Coeff0", &Coeff0);
        chain->SetBranchAddress("Coeff1", &Coeff1);
        chain->SetBranchAddress("Coeff2", &Coeff2);
        chain->SetBranchAddress("Coeff3", &Coeff3);
        chain->SetBranchAddress("Coeff4", &Coeff4);
        chain->SetBranchAddress("Coeff5", &Coeff5);
        chain->SetBranchAddress("Chi2", &Chi2);

        LogInfo << "TtRecoFile has " << chain->GetEntries() << " entries\n";

        return true;
    }

    bool find(const TimeStamp& ts_) {
        if (first_search) {
            return find_first(ts_);
        }

        Long64_t nentries = chain->GetEntries();
        if (nentries == 0l) return false;

        TimeStamp lower_bound = ts_ - TimeStamp{0, 1000};
        TimeStamp upper_bound = ts_;
        upper_bound.Add(TimeStamp{0, 1000});

        for (; cur_idx < nentries; ++cur_idx) {
            chain->GetEntry(cur_idx);
            
            TimeStamp cur_ts{start_TS->GetTimeSpec()};
            if (cur_ts < lower_bound) continue;
            else if (upper_bound < cur_ts) break;

            return true;
        }
        return false;
    }

    bool find_first(const TimeStamp& ts_) {
        Long64_t nentries = chain->GetEntries();
        if (nentries == 0l) return false;

        TimeStamp lower_bound = ts_ - TimeStamp{0, 1000};
        TimeStamp upper_bound = ts_;
        upper_bound.Add(TimeStamp{0, 1000});

        Long64_t left = 0l;
        Long64_t right = nentries - 1l;
        Long64_t result = -1l;

        while (left <= right) {
            Long64_t mid = left + (right - left) / 2;
            chain->GetEntry(mid);
            TimeStamp cur_ts{start_TS->GetTimeSpec()};
            
            if (cur_ts < lower_bound) {
                left = mid + 1;
            }
            else {
                result = mid;
                right = mid - 1;
            }
        
        }
        first_search = false;
        if (result == -1l) return false;
        
        cur_idx = result;
        chain->GetEntry(cur_idx);
        TimeStamp cur_ts{start_TS->GetTimeSpec()};

        return (lower_bound <= cur_ts) && (cur_ts <= upper_bound);
    }

};

struct TrackSaver {

    std::string filename;
    std::string treename = "muons";
    TFile* file = nullptr;
    TTree* tree = nullptr;

    int run_id = 0;
    time_t sec = 0l;
    int nsec = 0;
    double totq_cd = 0.0;
    double totq_wp = 0.0;

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
        tree = new TTree(treename.c_str(), treename.c_str());
        if (!tree) {
            LogError << "Cannot create track TTree\n";
            return false;
        }
        
        tree->Branch("run_id", &run_id);
        tree->Branch("sec", &sec);
        tree->Branch("nsec", &nsec);
        tree->Branch("totq_cd", &totq_cd);
        tree->Branch("totq_wp", &totq_wp);
        
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
        totq_cd = 0.0;
        totq_wp = 0.0;
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

    void fill(JM::EvtNavigator* nav) {
        if (!tree) return;
        std::shared_ptr<Event> evt = EventCache::load(nav);
        if (evt->tracks.empty()) return;
        reset();
        run_id = evt->run_id;
        sec = evt->ts.GetSec();
        nsec = evt->ts.GetNanoSec();
        totq_cd = evt->totq_cd;
        totq_wp = evt->totq_wp;
        for (const track& t : evt->tracks) {
            method.push_back(t.method);
            det.push_back(static_cast<unsigned char>(t.det));
            quality.push_back(t.quality);
            iposx.push_back(t.ipos.x);
            iposy.push_back(t.ipos.y);
            iposz.push_back(t.ipos.z);
            fposx.push_back(t.fpos.x);
            fposy.push_back(t.fpos.y);
            fposz.push_back(t.fpos.z);
        }
        if (method.empty()) return;
        tree->Fill();
    }

    bool save() {
        if (!file || !tree) return true;
        file->cd();
        tree->Write();
        file->Close();
        return true;
    }

};

struct FeatureSaver {

    std::string filename;
    std::string treename = "Features";
    TFile* file = nullptr;
    TTree* tree = nullptr;

    int run_id;
    time_t sec;
    int nsec;

    std::vector<double> iposx;
    std::vector<double> iposy;
    std::vector<double> iposz;
    std::vector<double> fposx;
    std::vector<double> fposy;
    std::vector<double> fposz;
    std::vector<double> chi2;
    std::vector<unsigned char> det; // 1 = CdWpTtChi2, 2 = CdClassify, 4 = TT

    std::vector<unsigned int> id;
    std::vector<double> fht;
    std::vector<double> totq;
    std::vector<double> q;
    std::vector<int> nhit;

    std::vector<double> pointx;
    std::vector<double> pointy;
    std::vector<double> pointz;

    bool init() {
        file = TFile::Open(filename.c_str(), "RECREATE");
        if (!file) {
            LogWarn << "Cannot open ROOT file " << filename << ". Skipping feature saving\n";
            return true;
        }
        tree = new TTree(treename.c_str(), treename.c_str());

        tree->Branch("run_id", &run_id);
        tree->Branch("sec", &sec);
        tree->Branch("nsec", &nsec);
        
        tree->Branch("iposx", &iposx);
        tree->Branch("iposy", &iposy);
        tree->Branch("iposz", &iposz);
        tree->Branch("fposx", &fposx);
        tree->Branch("fposy", &fposy);
        tree->Branch("fposz", &fposz);
        tree->Branch("chi2", &chi2);
        tree->Branch("det", &det);

        tree->Branch("id", &id);
        tree->Branch("fht", &fht);
        tree->Branch("totq", &totq);
        tree->Branch("q", &q);
        tree->Branch("nhit", &nhit);

        tree->Branch("pointx", &pointx);
        tree->Branch("pointy", &pointy);
        tree->Branch("pointz", &pointz);

        return true;  
    }
    
    void reset() {
        run_id = 0;
        sec = 0l;
        nsec = 0;

        iposx.clear();
        iposy.clear();
        iposz.clear();
        fposx.clear();
        fposy.clear();
        fposz.clear();
        chi2.clear();
        det.clear();

        id.clear();
        fht.clear();
        totq.clear();
        q.clear();
        nhit.clear();

        pointx.clear();
        pointy.clear();
        pointz.clear();
    }

    void fill() {
        if (tree && !chi2.empty() && !id.empty()) tree->Fill();
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

    double m_cd_muon_totq_thold = 10000.0;
    double m_cd_only_muon_totq_thold = 30000.0;
    double m_wp_muon_totq_thold = 600.0;
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
    FeatureSaver m_featureSaver;

	bool initBufSvc();
    bool initRecTool();
    bool initLoader();
    bool initAnalyses();

    calibration_context getCalibrationContext(const std::list<JM::CalibPmtChannel*>& clb_list);
    DetectorType getDetectorType(JM::EvtNavigator* nav);
    int getTtLayerId(double z);
    void addTtToTrack(std::vector<track>& tracks, const TimeStamp& curts);
    void addFeature(const std::vector<track>& tracks, const TimeStamp& curts, int run_id);

    void addTrack(RecTrks& rec_tracks, const std::string& method, double totq_cd, double totq_wp, const TimeStamp& ts, const track::loc& det, std::vector<track>& tracks);
    void addTrack(JM::CdTrackRecHeader* cdt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks);
    void addTrack(JM::WpRecHeader* wpt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks);
    void addTrack(JM::TtRecHeader* ttt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks);
    void addVertex(JM::OecHeader* oec_hdr, const std::string& method, const TimeStamp& ts, const calibration_context& calib, std::vector<vertex>& vertices);
    void addVertex(JM::CdVertexRecHeader* cdv_hdr, const std::string& method, const TimeStamp& ts, const calibration_context& calib, std::vector<vertex>& vertices);

};

#endif // ANALYSISGROUPC_ANALYSISGROUPC_HPP_