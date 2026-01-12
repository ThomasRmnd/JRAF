#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>

#include <TFile.h>
#include <TTree.h>

struct InputVars {

    int run_id;

    double posx_p;
    double posy_p;
    double posz_p;
    double e_p;
    time_t sec_p;
    int nsec_p;

    double totq_p;
    double meanq_p;
    double stdq_p;
    double minq_p;
    double maxq_p;
    std::size_t npmt_p;
    double meant_p;
    double stdt_p;

    double posx_d;
    double posy_d;
    double posz_d;
    double e_d;
    time_t sec_d;
    int nsec_d;
    
    double totq_d;
    double meanq_d;
    double stdq_d;
    double minq_d;
    double maxq_d;
    std::size_t npmt_d;
    double meant_d;
    double stdt_d;

    std::vector<double> *posx_n = nullptr;
    std::vector<double> *posy_n = nullptr;
    std::vector<double> *posz_n = nullptr;
    std::vector<double> *e_n = nullptr;
    std::vector<time_t> *sec_n = nullptr;
    std::vector<int> *nsec_n = nullptr;

    std::vector<double> *totq_n = nullptr;
    std::vector<double> *meanq_n = nullptr;
    std::vector<double> *stdq_n = nullptr;
    std::vector<double> *minq_n = nullptr;
    std::vector<double> *maxq_n = nullptr;
    std::vector<std::size_t> *npmt_n = nullptr;
    std::vector<double> *meant_n = nullptr;
    std::vector<double> *stdt_n = nullptr;

    std::vector<double> *posx_mult = nullptr;
    std::vector<double> *posy_mult = nullptr;
    std::vector<double> *posz_mult = nullptr;
    std::vector<double> *e_mult = nullptr;
    std::vector<time_t> *sec_mult = nullptr;
    std::vector<int> *nsec_mult = nullptr;
    std::vector<int> *mult_type = nullptr; // 0 = before prompt, 1 = between, 2 = after delayed

    std::vector<double> *totq_mult = nullptr;
    std::vector<double> *meanq_mult = nullptr;
    std::vector<double> *stdq_mult = nullptr;
    std::vector<double> *minq_mult = nullptr;
    std::vector<double> *maxq_mult = nullptr;
    std::vector<std::size_t> *npmt_mult = nullptr;
    std::vector<double> *meant_mult = nullptr;
    std::vector<double> *stdt_mult = nullptr;

    std::vector<std::string> *method_mu = nullptr;
    std::vector<int> *loc_mu = nullptr;
    std::vector<double> *posx_mu = nullptr;
    std::vector<double> *posy_mu = nullptr;
    std::vector<double> *posz_mu = nullptr;
    std::vector<double> *dirx_mu = nullptr;
    std::vector<double> *diry_mu = nullptr;
    std::vector<double> *dirz_mu = nullptr;
    std::vector<double> *totq_mu = nullptr;
    std::vector<time_t> *sec_mu = nullptr;
    std::vector<int> *nsec_mu = nullptr;
    std::vector<double> *quality_mu = nullptr;

};

struct OuputVars {

    int run_id;

    double posx_p;
    double posy_p;
    double posz_p;
    double e_p;
    time_t sec_p;
    int nsec_p;

    double totq_p;
    double meanq_p;
    double stdq_p;
    double minq_p;
    double maxq_p;
    std::size_t npmt_p;
    double meant_p;
    double stdt_p;

    double posx_d;
    double posy_d;
    double posz_d;
    double e_d;
    time_t sec_d;
    int nsec_d;
    
    double totq_d;
    double meanq_d;
    double stdq_d;
    double minq_d;
    double maxq_d;
    std::size_t npmt_d;
    double meant_d;
    double stdt_d;

    std::vector<double> posx_n;
    std::vector<double> posy_n;
    std::vector<double> posz_n;
    std::vector<double> e_n;
    std::vector<time_t> sec_n;
    std::vector<int> nsec_n;

    std::vector<double> totq_n;
    std::vector<double> meanq_n;
    std::vector<double> stdq_n;
    std::vector<double> minq_n;
    std::vector<double> maxq_n;
    std::vector<std::size_t> npmt_n;
    std::vector<double> meant_n;
    std::vector<double> stdt_n;

    std::vector<double> posx_mult;
    std::vector<double> posy_mult;
    std::vector<double> posz_mult;
    std::vector<double> e_mult;
    std::vector<time_t> sec_mult;
    std::vector<int> nsec_mult;
    std::vector<int> mult_type; // 0 = before prompt, 1 = between, 2 = after delayed

    std::vector<double> totq_mult;
    std::vector<double> meanq_mult;
    std::vector<double> stdq_mult;
    std::vector<double> minq_mult;
    std::vector<double> maxq_mult;
    std::vector<std::size_t> npmt_mult;
    std::vector<double> meant_mult;
    std::vector<double> stdt_mult;

    std::vector<std::string> method_mu;
    std::vector<int> loc_mu;
    std::vector<double> posx_mu;
    std::vector<double> posy_mu;
    std::vector<double> posz_mu;
    std::vector<double> dirx_mu;
    std::vector<double> diry_mu;
    std::vector<double> dirz_mu;
    std::vector<double> totq_mu;
    std::vector<time_t> sec_mu;
    std::vector<int> nsec_mu;
    std::vector<double> quality_mu;

};

void tree_fixer(const std::string& ifilename, const std::string& ofilename) {
    TFile* fin = TFile::Open(ifilename.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: Could not open input file: " << ifilename << '\n';
        return;
    }
    TFile* fout = TFile::Open(ofilename.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: Could not open output file: " << ofilename << '\n';
        return;
    }

    std::vector<std::string> treenames = {"IBDAnalysis__OMILREC", "IBDAnalysis__MixedPhase", "IBDAnalysis__OMILREC_JVtx"};

    for (const std::string& tname : treenames) {
        TTree* tin = dynamic_cast<TTree*>(fin->Get(tname.c_str()));
        if (!tin) {
            std::cerr << "ERROR: Could not find tree " << tname << " in file " << ifilename << '\n';
            fin->Close();
            return;
        }

        TTree* tout = tin->CloneTree(0);
        tout->SetName(tname.c_str());

        InputVars varin;
        tin->SetBranchAddress("posx_p", &varin.posx_p);
        tin->SetBranchAddress("posy_p", &varin.posy_p);
        tin->SetBranchAddress("posz_p", &varin.posz_p);
        tin->SetBranchAddress("e_p", &varin.e_p);
        tin->SetBranchAddress("sec_p", &varin.sec_p);
        tin->SetBranchAddress("nsec_p", &varin.nsec_p);
        tin->SetBranchAddress("totq_p", &varin.totq_p);
        tin->SetBranchAddress("meanq_p", &varin.meanq_p);
        tin->SetBranchAddress("stdq_p", &varin.stdq_p);
        tin->SetBranchAddress("minq_p", &varin.minq_p);
        tin->SetBranchAddress("maxq_p", &varin.maxq_p);
        tin->SetBranchAddress("npmt_p", &varin.npmt_p);
        tin->SetBranchAddress("meant_p", &varin.meant_p);
        tin->SetBranchAddress("stdt_p", &varin.stdt_p);
        tin->SetBranchAddress("posx_d", &varin.posx_d);
        tin->SetBranchAddress("posy_d", &varin.posy_d);
        tin->SetBranchAddress("posz_d", &varin.posz_d);
        tin->SetBranchAddress("e_d", &varin.e_d);
        tin->SetBranchAddress("sec_d", &varin.sec_d);
        tin->SetBranchAddress("nsec_d", &varin.nsec_d);
        tin->SetBranchAddress("totq_d", &varin.totq_d);
        tin->SetBranchAddress("meanq_d", &varin.meanq_d);
        tin->SetBranchAddress("stdq_d", &varin.stdq_d);
        tin->SetBranchAddress("minq_d", &varin.minq_d);
        tin->SetBranchAddress("maxq_d", &varin.maxq_d);
        tin->SetBranchAddress("npmt_d", &varin.npmt_d);
        tin->SetBranchAddress("meant_d", &varin.meant_d);
        tin->SetBranchAddress("stdt_d", &varin.stdt_d);
        tin->SetBranchAddress("posx_n", &varin.posx_n);
        tin->SetBranchAddress("posy_n", &varin.posy_n);
        tin->SetBranchAddress("posz_n", &varin.posz_n);
        tin->SetBranchAddress("e_n", &varin.e_n);
        tin->SetBranchAddress("sec_n", &varin.sec_n);
        tin->SetBranchAddress("nsec_n", &varin.nsec_n);
        tin->SetBranchAddress("totq_n", &varin.totq_n);
        tin->SetBranchAddress("meanq_n", &varin.meanq_n);
        tin->SetBranchAddress("stdq_n", &varin.stdq_n);
        tin->SetBranchAddress("minq_n", &varin.minq_n);
        tin->SetBranchAddress("maxq_n", &varin.maxq_n);
        tin->SetBranchAddress("npmt_n", &varin.npmt_n);
        tin->SetBranchAddress("meant_n", &varin.meant_n);
        tin->SetBranchAddress("stdt_n", &varin.stdt_n);
        tin->SetBranchAddress("posx_mult", &varin.posx_mult);
        tin->SetBranchAddress("posy_mult", &varin.posy_mult);
        tin->SetBranchAddress("posz_mult", &varin.posz_mult);
        tin->SetBranchAddress("e_mult", &varin.e_mult);
        tin->SetBranchAddress("sec_mult", &varin.sec_mult);
        tin->SetBranchAddress("nsec_mult", &varin.nsec_mult);
        tin->SetBranchAddress("mult_type", &varin.mult_type);
        tin->SetBranchAddress("totq_mult", &varin.totq_mult);
        tin->SetBranchAddress("meanq_mult", &varin.meanq_mult);
        tin->SetBranchAddress("stdq_mult", &varin.stdq_mult);
        tin->SetBranchAddress("minq_mult", &varin.minq_mult);
        tin->SetBranchAddress("maxq_mult", &varin.maxq_mult);
        tin->SetBranchAddress("npmt_mult", &varin.npmt_mult);
        tin->SetBranchAddress("meant_mult", &varin.meant_mult);
        tin->SetBranchAddress("stdt_mult", &varin.stdt_mult);
        tin->SetBranchAddress("method_mu", &varin.method_mu);
        tin->SetBranchAddress("loc_mu", &varin.loc_mu);
        tin->SetBranchAddress("posx_mu", &varin.posx_mu);
        tin->SetBranchAddress("posy_mu", &varin.posy_mu);
        tin->SetBranchAddress("posz_mu", &varin.posz_mu);
        tin->SetBranchAddress("dirx_mu", &varin.dirx_mu);
        tin->SetBranchAddress("diry_mu", &varin.diry_mu);
        tin->SetBranchAddress("dirz_mu", &varin.dirz_mu);
        tin->SetBranchAddress("totq_mu", &varin.totq_mu);
        tin->SetBranchAddress("sec_mu", &varin.sec_mu);
        tin->SetBranchAddress("nsec_mu", &varin.nsec_mu);
        tin->SetBranchAddress("quality_mu", &varin.quality_mu);

        OuputVars varout;
        tout->Branch("posx_p", &varout.posx_p);
        tout->Branch("posy_p", &varout.posy_p);
        tout->Branch("posz_p", &varout.posz_p);
        tout->Branch("e_p", &varout.e_p);
        tout->Branch("sec_p", &varout.sec_p);
        tout->Branch("nsec_p", &varout.nsec_p);
        tout->Branch("totq_p", &varout.totq_p);
        tout->Branch("meanq_p", &varout.meanq_p);
        tout->Branch("stdq_p", &varout.stdq_p);
        tout->Branch("minq_p", &varout.minq_p);
        tout->Branch("maxq_p", &varout.maxq_p);
        tout->Branch("npmt_p", &varout.npmt_p);
        tout->Branch("meant_p", &varout.meant_p);
        tout->Branch("stdt_p", &varout.stdt_p);
        tout->Branch("posx_d", &varout.posx_d);
        tout->Branch("posy_d", &varout.posy_d);
        tout->Branch("posz_d", &varout.posz_d);
        tout->Branch("e_d", &varout.e_d);
        tout->Branch("sec_d", &varout.sec_d);
        tout->Branch("nsec_d", &varout.nsec_d);
        tout->Branch("totq_d", &varout.totq_d);
        tout->Branch("meanq_d", &varout.meanq_d);
        tout->Branch("stdq_d", &varout.stdq_d);
        tout->Branch("minq_d", &varout.minq_d);
        tout->Branch("maxq_d", &varout.maxq_d);
        tout->Branch("npmt_d", &varout.npmt_d);
        tout->Branch("meant_d", &varout.meant_d);
        tout->Branch("stdt_d", &varout.stdt_d);
        tout->Branch("posx_n", &varout.posx_n);
        tout->Branch("posy_n", &varout.posy_n);
        tout->Branch("posz_n", &varout.posz_n);
        tout->Branch("e_n", &varout.e_n);
        tout->Branch("sec_n", &varout.sec_n);
        tout->Branch("nsec_n", &varout.nsec_n);
        tout->Branch("totq_n", &varout.totq_n);
        tout->Branch("meanq_n", &varout.meanq_n);
        tout->Branch("stdq_n", &varout.stdq_n);
        tout->Branch("minq_n", &varout.minq_n);
        tout->Branch("maxq_n", &varout.maxq_n);
        tout->Branch("npmt_n", &varout.npmt_n);
        tout->Branch("meant_n", &varout.meant_n);
        tout->Branch("stdt_n", &varout.stdt_n);
        tout->Branch("posx_mult", &varout.posx_mult);
        tout->Branch("posy_mult", &varout.posy_mult);
        tout->Branch("posz_mult", &varout.posz_mult);
        tout->Branch("e_mult", &varout.e_mult);
        tout->Branch("sec_mult", &varout.sec_mult);
        tout->Branch("nsec_mult", &varout.nsec_mult);
        tout->Branch("mult_type", &varout.mult_type);
        tout->Branch("totq_mult", &varout.totq_mult);
        tout->Branch("meanq_mult", &varout.meanq_mult);
        tout->Branch("stdq_mult", &varout.stdq_mult);
        tout->Branch("minq_mult", &varout.minq_mult);
        tout->Branch("maxq_mult", &varout.maxq_mult);
        tout->Branch("npmt_mult", &varout.npmt_mult);
        tout->Branch("meant_mult", &varout.meant_mult);
        tout->Branch("stdt_mult", &varout.stdt_mult);
        tout->Branch("method_mu", &varout.method_mu);
        tout->Branch("loc_mu", &varout.loc_mu);
        tout->Branch("posx_mu", &varout.posx_mu);
        tout->Branch("posy_mu", &varout.posy_mu);
        tout->Branch("posz_mu", &varout.posz_mu);
        tout->Branch("dirx_mu", &varout.dirx_mu);
        tout->Branch("diry_mu", &varout.diry_mu);
        tout->Branch("dirz_mu", &varout.dirz_mu);
        tout->Branch("totq_mu", &varout.totq_mu);
        tout->Branch("sec_mu", &varout.sec_mu);
        tout->Branch("nsec_mu", &varout.nsec_mu);
        tout->Branch("quality_mu", &varout.quality_mu);

        Long64_t nentries = tin->GetEntries();
    
        std::cout << "Starting fix for " << nentries << " events in " << tname << "...\n";

        for (Long64_t i = 0; i < nentries; ++i) {
            tin->GetEntry(i);

            varout.run_id = varin.run_id;
            varout.posx_p = varin.posx_p;
            varout.posy_p = varin.posy_p;
            varout.posz_p = varin.posz_p;
            varout.e_p = varin.e_p;
            varout.sec_p = varin.sec_p;
            varout.nsec_p = varin.nsec_p;
            varout.totq_p = varin.totq_p;
            varout.meanq_p = varin.meanq_p;
            varout.stdq_p = varin.stdq_p;
            varout.minq_p = varin.minq_p;
            varout.maxq_p = varin.maxq_p;
            varout.npmt_p = varin.npmt_p;
            varout.meant_p = varin.meant_p;
            varout.stdt_p = varin.stdt_p;
            varout.posx_d = varin.posx_d;
            varout.posy_d = varin.posy_d;
            varout.posz_d = varin.posz_d;
            varout.e_d = varin.e_d;
            varout.sec_d = varin.sec_d;
            varout.nsec_d = varin.nsec_d;
            varout.totq_d = varin.totq_d;
            varout.meanq_d = varin.meanq_d;
            varout.stdq_d = varin.stdq_d;
            varout.minq_d = varin.minq_d;
            varout.maxq_d = varin.maxq_d;
            varout.npmt_d = varin.npmt_d;
            varout.meant_d = varin.meant_d;
            varout.stdt_d = varin.stdt_d;

            varout.method_mu = *varin.method_mu;
            varout.loc_mu = *varin.loc_mu;
            varout.posx_mu = *varin.posx_mu;
            varout.posy_mu = *varin.posy_mu;
            varout.posz_mu = *varin.posz_mu;
            varout.dirx_mu = *varin.dirx_mu;
            varout.diry_mu = *varin.diry_mu;
            varout.dirz_mu = *varin.dirz_mu;
            varout.totq_mu = *varin.totq_mu;
            varout.sec_mu = *varin.sec_mu;
            varout.nsec_mu = *varin.nsec_mu;
            varout.quality_mu = *varin.quality_mu;

            varout.posx_n = *varin.posx_n;
            varout.posy_n = *varin.posy_n;
            varout.posz_n = *varin.posz_n;
            varout.e_n = *varin.e_n;
            varout.sec_n = *varin.sec_n;
            varout.nsec_n = *varin.nsec_n;

            varout.totq_n = *varin.totq_n;
            varout.meanq_n.assign(varin.meanq_n->end() - varin.e_n->size(), varin.meanq_n->end());
            varout.stdq_n.assign(varin.stdq_n->end() - varin.e_n->size(), varin.stdq_n->end());
            varout.minq_n.assign(varin.minq_n->end() - varin.e_n->size(), varin.minq_n->end());
            varout.maxq_n.assign(varin.maxq_n->end() - varin.e_n->size(), varin.maxq_n->end());
            varout.npmt_n.assign(varin.npmt_n->end() - varin.e_n->size(), varin.npmt_n->end());
            varout.meant_n.assign(varin.meant_n->end() - varin.e_n->size(), varin.meant_n->end());
            varout.stdt_n.assign(varin.stdt_n->end() - varin.e_n->size(), varin.stdt_n->end());

            varout.posx_mult = *varin.posx_mult;
            varout.posy_mult = *varin.posy_mult;
            varout.posz_mult = *varin.posz_mult;
            varout.e_mult = *varin.e_mult;
            varout.sec_mult = *varin.sec_mult;
            varout.nsec_mult = *varin.nsec_mult;
            varout.mult_type = *varin.mult_type;

            varout.totq_mult = *varin.totq_mult;
            varout.meanq_mult.assign(varin.meanq_mult->end() - varin.e_mult->size(), varin.meanq_mult->end());
            varout.stdq_mult.assign(varin.stdq_mult->end() - varin.e_mult->size(), varin.stdq_mult->end());
            varout.minq_mult.assign(varin.minq_mult->end() - varin.e_mult->size(), varin.minq_mult->end());
            varout.maxq_mult.assign(varin.maxq_mult->end() - varin.e_mult->size(), varin.maxq_mult->end());
            varout.npmt_mult.assign(varin.npmt_mult->end() - varin.e_mult->size(), varin.npmt_mult->end());
            varout.meant_mult.assign(varin.meant_mult->end() - varin.e_mult->size(), varin.meant_mult->end());
            varout.stdt_mult.assign(varin.stdt_mult->end() - varin.e_mult->size(), varin.stdt_mult->end());

            tout->Fill();

            if (i % 1000 == 0) {
                std::cout << "Processing (" << tname << ") event " << i << "...\n";
            }
        }

        std::cout << "Fixing complete. Writing and closing files\n";

        fout->cd();
        tout->Write();
        delete tout;
    }
    fin->Close();
    fout->Close();
}