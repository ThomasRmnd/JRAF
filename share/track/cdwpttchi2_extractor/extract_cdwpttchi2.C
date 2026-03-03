#include <iostream>

#include <TFile.h>
#include <TTree.h>

struct IBranches {
    int run_id;
    time_t sec;
    int nsec;
    double totq_cd;
    double totq_wp;
    std::vector<std::string> *method = nullptr;
    std::vector<unsigned char> *det = nullptr;
    std::vector<double> *quality = nullptr;
    std::vector<double> *iposx = nullptr, *iposy = nullptr, *iposz = nullptr; 
    std::vector<double> *fposx = nullptr, *fposy = nullptr, *fposz = nullptr;
};

struct OBranches {
    int run_id;
    time_t sec;
    int nsec;
    double totq_cd;
    double totq_wp;
    double chi2;
    double iposx, iposy, iposz;
    double fposx, fposy, fposz;
    unsigned int ntracks_cdclassify;
    unsigned int ntracks_wpclassify;
    unsigned int nstoppings_cdclassify;
    unsigned int nstoppings_wpclassify;
};

int extract_cdwpttchi2(const char* ipath, const char* opath) {
    TFile* ifile = TFile::Open(ipath, "READ");
    if (!ifile) {
        std::cerr << "Error: Unable to open file " << ipath << '\n';
        return 1;
    }
    TTree* itree = ifile->Get<TTree>("muons");
    if (!itree) {
        std::cerr << "Error: Cannot retrieve tree muons in file " << ipath << '\n';
        return 1;
    }
    TFile* ofile = TFile::Open(opath, "RECREATE");
    if (!ofile) {
        std::cerr << "Error: Unable to open file " << opath << '\n';
        return 1;
    }
    TTree* otree = new TTree("muons", "muons");
    if (!otree) {
        std::cerr << "Error: Unable to create output tree\n";
        return 1;
    }

    IBranches ib;
    itree->SetBranchAddress("run_id", &ib.run_id);
    itree->SetBranchAddress("sec", &ib.sec);
    itree->SetBranchAddress("nsec", &ib.nsec);
    itree->SetBranchAddress("totq_cd", &ib.totq_cd);
    itree->SetBranchAddress("totq_wp", &ib.totq_wp);
    itree->SetBranchAddress("method", &ib.method);
    itree->SetBranchAddress("det", &ib.det);
    itree->SetBranchAddress("quality", &ib.quality);
    itree->SetBranchAddress("iposx", &ib.iposx);
    itree->SetBranchAddress("iposy", &ib.iposy);
    itree->SetBranchAddress("iposz", &ib.iposz);
    itree->SetBranchAddress("fposx", &ib.fposx);
    itree->SetBranchAddress("fposy", &ib.fposy);
    itree->SetBranchAddress("fposz", &ib.fposz);

    OBranches ob;
    otree->Branch("run_id", &ob.run_id);
    otree->Branch("sec", &ob.sec);
    otree->Branch("nsec", &ob.nsec);
    otree->Branch("totq_cd", &ob.totq_cd);
    otree->Branch("totq_wp", &ob.totq_wp);
    otree->Branch("chi2", &ob.chi2);
    otree->Branch("iposx", &ob.iposx);
    otree->Branch("iposy", &ob.iposy);
    otree->Branch("iposz", &ob.iposz);
    otree->Branch("fposx", &ob.fposx);
    otree->Branch("fposy", &ob.fposy);
    otree->Branch("fposz", &ob.fposz);
    otree->Branch("ntracks_cdclassify", &ob.ntracks_cdclassify);
    otree->Branch("ntracks_wpclassify", &ob.ntracks_wpclassify);
    otree->Branch("nstoppings_cdclassify", &ob.nstoppings_cdclassify);
    otree->Branch("nstoppings_wpclassify", &ob.nstoppings_wpclassify);

    long nentries = itree->GetEntries();
    std::cout << "Info: Found " << nentries << " entries in input file\n";

    for (long k = 0l; k < nentries; ++k) {
        itree->GetEntry(k);
        bool found = false;
        ob.run_id = ib.run_id;
        ob.sec = ib.sec;
        ob.nsec = ib.nsec;
        ob.totq_cd = ib.totq_cd;
        ob.totq_wp = ib.totq_wp;
        unsigned int ntracks_cdclassify = 0u;
        unsigned int ntracks_wpclassify = 0u;
        unsigned int nstoppings_cdclassify = 0u;
        unsigned int nstoppings_wpclassify = 0u;
        for (std::size_t i = 0ul; i < ib.method->size(); ++i) {
            if ((*ib.method)[i] == "CdClassify") {
                ++ntracks_cdclassify;
                TVector3 fpos((*ib.fposx)[i], (*ib.fposy)[i], (*ib.fposz)[i]);
                if (fpos.Mag() > 40000.0) {
                    ++nstoppings_cdclassify;
                }
            } 
            else if ((*ib.method)[i] == "WpBasic") {
                ++ntracks_wpclassify;
                TVector3 fpos((*ib.fposx)[i], (*ib.fposy)[i], (*ib.fposz)[i]);
                if (fpos.Mag() > 40000.0) {
                    ++nstoppings_wpclassify;
                }
            }
            else if ((*ib.method)[i] == "CdWpTtChi2") {
                ob.chi2 = (*ib.quality)[i];
                ob.iposx = (*ib.iposx)[i];
                ob.iposy = (*ib.iposy)[i];
                ob.iposz = (*ib.iposz)[i];
                ob.fposx = (*ib.fposx)[i];
                ob.fposy = (*ib.fposy)[i];
                ob.fposz = (*ib.fposz)[i];
                if (std::isnan(ob.chi2) || 
                    std::isnan(ob.iposx) || std::isnan(ob.iposy) || std::isnan(ob.iposz) || 
                    std::isnan(ob.fposx) || std::isnan(ob.fposy) || std::isnan(ob.fposz)
                ) {
                    continue;
                }
                found = true;
            }
        }
        ob.ntracks_cdclassify = ntracks_cdclassify;
        ob.ntracks_wpclassify = ntracks_wpclassify;
        ob.nstoppings_cdclassify = nstoppings_cdclassify;
        ob.nstoppings_wpclassify = nstoppings_wpclassify;
        if (found) {
            otree->Fill();
        }
    }

    ofile->Write();
    ofile->Close();
    ifile->Close();

    return 0;
}