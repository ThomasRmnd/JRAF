#ifndef JRAF_UTILS_TTRECOFILE_HPP_
#define JRAF_UTILS_TTRECOFILE_HPP_

#include <string>

#include <TChain.h>
#include <TTimeStamp.h>

#include "SniperKernel/SniperLog.h"

#include "Context/TimeStamp.h"

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

#endif // JRAF_UTILS_TTRECOFILE_HPP_