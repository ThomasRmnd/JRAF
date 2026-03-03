#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>

TTimeStamp operator+(const TTimeStamp& lhs, const TTimeStamp& rhs) {
    return TTimeStamp(lhs.GetSec() + rhs.GetSec(), lhs.GetNanoSec() + rhs.GetNanoSec());
}

TTimeStamp operator-(const TTimeStamp& lhs, const TTimeStamp& rhs) {
    return TTimeStamp(lhs.GetSec() - rhs.GetSec(), lhs.GetNanoSec() - rhs.GetNanoSec());
}

struct vertex {
    TVector3 pos;
    TTimeStamp ts;
    double e;
};

bool operator<(const vertex& lhs, const vertex& rhs) {
    return lhs.ts < rhs.ts;
}

struct ibd {
    int run_id;
    vertex prompt;
    vertex delayed;
};

bool operator<(const ibd& lhs, const ibd& rhs) {
    return lhs.prompt < rhs.prompt;
}

struct track {
    std::string method;
    TVector3 ipos;
    TVector3 fpos;
    double quality;
    unsigned char det;
};

struct muon {
    int run_id;
    TTimeStamp ts;
    double totq_cd;
    double totq_wp;
    std::vector<track> tracks;
};

bool operator<(const muon& lhs, const muon& rhs) {
    return lhs.ts < rhs.ts;
}

int estimate_cosmo_rate(const char* datapath, const char* muonpath) {
    TFile* file_data = TFile::Open(datapath, "READ");
    if (!file_data) {
        std::cerr << "Cannot open file " << datapath << '\n';
        return 1;
    }
    TTree* tree_data = (TTree*)file_data->Get("events");
    if (!tree_data) {
        std::cerr << "Cannot get tree 'events' from file " << datapath << '\n';
        return 1;
    }
    TChain* chain_muon = new TChain("muons");
    chain_muon->Add(muonpath);
    if (chain_muon->GetEntries() == 0) {
        std::cerr << "Cannot open file " << muonpath << '\n';
        return 1;
    }

    int run_id;
    double posx_p, posy_p, posz_p;
    time_t sec_p;
    int nsec_p;
    double e_p;
    double posx_d, posy_d, posz_d;
    time_t sec_d;
    int nsec_d;
    double e_d;

    tree_data->SetBranchAddress("run_id", &run_id);
    tree_data->SetBranchAddress("posx_p", &posx_p);
    tree_data->SetBranchAddress("posy_p", &posy_p);
    tree_data->SetBranchAddress("posz_p", &posz_p);
    tree_data->SetBranchAddress("sec_p", &sec_p);
    tree_data->SetBranchAddress("nsec_p", &nsec_p);
    tree_data->SetBranchAddress("e_p", &e_p);
    tree_data->SetBranchAddress("posx_d", &posx_d);
    tree_data->SetBranchAddress("posy_d", &posy_d);
    tree_data->SetBranchAddress("posz_d", &posz_d);
    tree_data->SetBranchAddress("sec_d", &sec_d);
    tree_data->SetBranchAddress("nsec_d", &nsec_d);
    tree_data->SetBranchAddress("e_d", &e_d);

    std::vector<ibd> ibds;
    for (int i = 0; i < tree_data->GetEntries(); ++i) {
        tree_data->GetEntry(i);
        ibds.push_back({
            run_id,
            {
                TVector3(posx_p, posy_p, posz_p),
                TTimeStamp(sec_p, nsec_p),
                e_p
            },
            {
                TVector3(posx_d, posy_d, posz_d),
                TTimeStamp(sec_d, nsec_d),
                e_d
            }
        });
    }
    std::sort(ibds.begin(), ibds.end());
    std::cout << "Info: " << ibds.size() << " IBDs event retrieved\n";

    time_t sec;
    int nsec;
    double totq_cd, totq_wp;
    std::vector<std::string>* method = nullptr;
    std::vector<unsigned char>* det = nullptr;
    std::vector<double>* quality = nullptr;
    std::vector<double> *iposx = nullptr, *iposy = nullptr, *iposz = nullptr;
    std::vector<double> *fposx = nullptr, *fposy = nullptr, *fposz = nullptr;

    chain_muon->SetBranchAddress("run_id", &run_id);
    chain_muon->SetBranchAddress("sec", &sec);
    chain_muon->SetBranchAddress("nsec", &nsec);
    chain_muon->SetBranchAddress("totq_cd", &totq_cd);
    chain_muon->SetBranchAddress("totq_wp", &totq_wp);
    chain_muon->SetBranchAddress("method", &method);
    chain_muon->SetBranchAddress("det", &det);
    chain_muon->SetBranchAddress("quality", &quality);
    chain_muon->SetBranchAddress("iposx", &iposx);
    chain_muon->SetBranchAddress("iposy", &iposy);
    chain_muon->SetBranchAddress("iposz", &iposz);
    chain_muon->SetBranchAddress("fposx", &fposx);
    chain_muon->SetBranchAddress("fposy", &fposy);
    chain_muon->SetBranchAddress("fposz", &fposz);

    std::vector<muon> muons;
    for (int i = 0; i < chain_muon->GetEntries(); ++i) {
        chain_muon->GetEntry(i);
        muon mu;
        mu.run_id = run_id;
        mu.ts = TTimeStamp(sec, nsec);
        mu.totq_cd = totq_cd;
        mu.totq_wp = totq_wp;
        for (int j = 0; j < method->size(); ++j) {
            mu.tracks.push_back({
                method->at(j),
                TVector3(iposx->at(j), iposy->at(j), iposz->at(j)),
                TVector3(fposx->at(j), fposy->at(j), fposz->at(j)),
                quality->at(j),
                det->at(j)
            });
        }
        muons.push_back(mu);
    }
    std::sort(muons.begin(), muons.end());
    std::cout << "Info: " << muons.size() << " muons event retrieved\n";

    std::size_t j = 0ul;
    for (const ibd& event : ibds) {
        std::vector<muon>::const_iterator it = std::upper_bound(muons.begin(), muons.end(), event.prompt.ts, [](const TTimeStamp& ts, const muon& mu) { return ts < mu.ts; });
        if (it == muons.begin()) {
            std::cout << "IBD has no prior muon\n";
            continue;
        }
        --it;
        TTimeStamp dt = event.prompt.ts - it->ts;
        std::cout << "dt = " << dt.AsDouble() << " s\n";
    }

    return 0;
}