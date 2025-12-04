#include <cmath>
#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TVector3.h>

inline TTimeStamp operator-(const TTimeStamp& lhs, const TTimeStamp& rhs) {
    return TTimeStamp(lhs.GetSec() - rhs.GetSec(), lhs.GetNanoSec() - rhs.GetNanoSec());
}

inline TTimeStamp operator+(const TTimeStamp& lhs, const TTimeStamp& rhs) {
    return TTimeStamp(lhs.GetSec() + rhs.GetSec(), lhs.GetNanoSec() + rhs.GetNanoSec());
}

void compare_with_TT(const char* cdwp_filename, const char* tt_filename, const char* output_filename) {
    TChain tt_chain("TT");
    tt_chain.Add(tt_filename);
    
    Int_t evtID, NTotPoints, NTracks;
    Int_t NPoints[20];
    TTimeStamp* start_TS = nullptr;
    Float_t PointX[20], PointY[20], PointZ[20];
    Double_t Coeff0[20], Coeff1[20], Coeff2[20], Coeff3[20], Coeff4[20], Coeff5[20];
    Double_t Chi2[20];

    tt_chain.SetBranchAddress("evtID", &evtID);
    tt_chain.SetBranchAddress("NTotPoints", &NTotPoints);
    tt_chain.SetBranchAddress("PointX", &PointX);
    tt_chain.SetBranchAddress("PointY", &PointY);
    tt_chain.SetBranchAddress("PointZ", &PointZ);
    tt_chain.SetBranchAddress("NTracks", &NTracks);
    tt_chain.SetBranchAddress("NPoints", NPoints);
    tt_chain.SetBranchAddress("start_TS", &start_TS);
    tt_chain.SetBranchAddress("Coeff0", &Coeff0);
    tt_chain.SetBranchAddress("Coeff1", &Coeff1);
    tt_chain.SetBranchAddress("Coeff2", &Coeff2);
    tt_chain.SetBranchAddress("Coeff3", &Coeff3);
    tt_chain.SetBranchAddress("Coeff4", &Coeff4);
    tt_chain.SetBranchAddress("Coeff5", &Coeff5);
    tt_chain.SetBranchAddress("Chi2", &Chi2);

    int cur_tt_idx = 0;

    TFile* cdwp_file = TFile::Open(cdwp_filename);
    if (!cdwp_file) {
        std::cerr << "Cannot open CD-WP muon reconstruction file\n";
    return;
    }

    TTree* cdwp_tree = cdwp_file->Get<TTree>("muons");
    if (!cdwp_tree) {
        std::cerr << "Cannot retrieve CD-WP muon reconstruction tree\n";
    return;
    }

    Int_t run_id;
    Long_t sec;
    Int_t nsec;
    std::vector<std::string>* method = nullptr;
    std::vector<unsigned char>* det = nullptr;
    std::vector<double>* quality = nullptr;
    std::vector<double> *iposx = nullptr, *iposy = nullptr, *iposz = nullptr, *fposx = nullptr, *fposy = nullptr, *fposz = nullptr;

    cdwp_tree->SetBranchAddress("run_id", &run_id);
    cdwp_tree->SetBranchAddress("sec", &sec);
    cdwp_tree->SetBranchAddress("nsec", &nsec);
    cdwp_tree->SetBranchAddress("method", &method);
    cdwp_tree->SetBranchAddress("det", &det);
    cdwp_tree->SetBranchAddress("quality", &quality);
    cdwp_tree->SetBranchAddress("iposx", &iposx);
    cdwp_tree->SetBranchAddress("iposy", &iposy);
    cdwp_tree->SetBranchAddress("iposz", &iposz);
    cdwp_tree->SetBranchAddress("fposx", &fposx);
    cdwp_tree->SetBranchAddress("fposy", &fposy);
    cdwp_tree->SetBranchAddress("fposz", &fposz);

    long tt_nentries = tt_chain.GetEntries();
    long cdwp_nentries = cdwp_tree->GetEntries();

    std::size_t nmatches = 0ul;

    std::cout << "CD-WP: " << cdwp_nentries << ", TT: " << tt_nentries << '\n';

    TFile* fout = TFile::Open(output_filename, "RECREATE");
    TTree* tree_out = new TTree("matches", "Matched CDWP and TT tracks");

    Int_t run_id_out;
    Long_t sec_out;
    Int_t nsec_out;
    std::vector<std::string> method_out;
    std::vector<unsigned char> det_out;
    std::vector<double> quality_out;
    std::vector<double> iposx_out, iposy_out, iposz_out, fposx_out, fposy_out, fposz_out;

    tree_out->Branch("run_id", &run_id_out);
    tree_out->Branch("sec", &sec_out);
    tree_out->Branch("nsec", &nsec_out);
    tree_out->Branch("method", &method_out);
    tree_out->Branch("det", &det_out);
    tree_out->Branch("quality", &quality_out);
    tree_out->Branch("iposx", &iposx_out);
    tree_out->Branch("iposy", &iposy_out);
    tree_out->Branch("iposz", &iposz_out);
    tree_out->Branch("fposx", &fposx_out);
    tree_out->Branch("fposy", &fposy_out);
    tree_out->Branch("fposz", &fposz_out);

    Int_t evtIDOut, NTotPointsOut, NTracksOut;
    Float_t PointXOut[20], PointYOut[20], PointZOut[20];
    Int_t NPointsOut[20];
    TTimeStamp start_TSOut;
    Double_t Coeff0Out[20], Coeff1Out[20], Coeff2Out[20], Coeff3Out[20], Coeff4Out[20], Coeff5Out[20];
    Double_t Chi2Out[20];

    tree_out->Branch("evtID", &evtIDOut);
    tree_out->Branch("NTotPoints", &NTotPointsOut, "NTotPoints/I");
    tree_out->Branch("PointX", &PointXOut, "PointX[20]/F");
    tree_out->Branch("PointY", &PointYOut, "PointY[20]/F");
    tree_out->Branch("PointZ", &PointZOut, "PointZ[20]/F");
    tree_out->Branch("NTracks", &NTracksOut, "NTracks/I");
    tree_out->Branch("start_TS", &start_TSOut);
    tree_out->Branch("NPoints", NPointsOut, "NPoints[20]/I");
    tree_out->Branch("Coeff0", Coeff0Out, "Coeff0[20]/D");
    tree_out->Branch("Coeff1", Coeff1Out, "Coeff1[20]/D");
    tree_out->Branch("Coeff2", Coeff2Out, "Coeff2[20]/D");
    tree_out->Branch("Coeff3", Coeff3Out, "Coeff3[20]/D");
    tree_out->Branch("Coeff4", Coeff4Out, "Coeff4[20]/D");
    tree_out->Branch("Coeff5", Coeff5Out, "Coeff5[20]/D");
    tree_out->Branch("Chi2", Chi2Out, "Chi2[20]/D");

    TVector3 cdwp_pos, tt_pos, cdwp_dir, tt_dir;

    for (int k = 0; k < cdwp_nentries; ++k) {
        cdwp_tree->GetEntry(k);

        TTimeStamp cdwp_ts{sec, nsec};
        TTimeStamp lower_bound = cdwp_ts - TTimeStamp{0, 1000};
        TTimeStamp upper_bound = cdwp_ts + TTimeStamp{0, 1000};

        cdwp_pos.SetXYZ(iposx->at(0), iposy->at(0), iposz->at(0));
        cdwp_dir.SetXYZ(fposx->at(0), fposy->at(0), fposz->at(0));
        cdwp_dir = (cdwp_dir - cdwp_pos).Unit();
        
        run_id_out = run_id;
        sec_out = sec;
        nsec_out = nsec;
        method_out = *method;
        det_out = *det;
        quality_out = *quality;
        iposx_out = *iposx;
        iposy_out = *iposy;
        iposz_out = *iposz;
        fposx_out = *fposx;
        fposy_out = *fposy;
        fposz_out = *fposz;

        evtIDOut = 0;
        NTotPointsOut = 0;
        NTracksOut = 0;
        start_TSOut = TTimeStamp{0, 0};

        for (; cur_tt_idx < tt_nentries; ++cur_tt_idx) {
            tt_chain.GetEntry(cur_tt_idx);

            TTimeStamp tt_ts{*start_TS};

            if (tt_ts < lower_bound) continue;
            if (upper_bound < tt_ts) break;

            ++nmatches;
            evtIDOut = evtID;
            NTotPointsOut = NTotPoints;
            NTracksOut = NTracks;
            start_TSOut = *start_TS;
            for (int i = 0; i < 20; ++i) {
                PointXOut[i] = PointX[i];
                PointYOut[i] = PointY[i];
                PointZOut[i] = PointZ[i];
                NPointsOut[i] = NPoints[i];
                Coeff0Out[i] = Coeff0[i];
                Coeff1Out[i] = Coeff1[i];
                Coeff2Out[i] = Coeff2[i];
                Coeff3Out[i] = Coeff3[i];
                Coeff4Out[i] = Coeff4[i];
                Coeff5Out[i] = Coeff5[i];
                Chi2Out[i] = Chi2[i];
            }

            std::cout << "Matching! TimeStamp: " << TTimeStamp{sec, nsec} << '\n';
            for (std::size_t k = 0ul; k < method->size(); ++k) {
                std::cout << "Method: " << (*method)[k] << ", Det: " << (*det)[k] << ", Quality: " << (*quality)[k] << ", ix: " << (*iposx)[k] << ", iy: " << (*iposy)[k] << ", iz: " << (*iposz)[k] << ", fx: " << (*fposx)[k] << ", fy: " << (*fposy)[k] << ", fz: " << (*fposz)[k] << '\n';
            }
            for (int k = 0; k < NTracks; ++k) {
                std::cout << "Chi2: " << Chi2[k] << ", x: " << Coeff0[k] << ", y: " << Coeff1[k] << ", z: " << Coeff2[k] << ", dx: " << Coeff3[k] << ", dy: " << Coeff4[k] << ", " << Coeff5[k] << '\n';
            }
        }

        tree_out->Fill();   
    }

    std::cout << "Total matches: " << nmatches << '\n';
    fout->cd();
    tree_out->Write();
    fout->Close();
}