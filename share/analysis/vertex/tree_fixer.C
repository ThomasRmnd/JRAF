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
    std::size_t nhit_p;
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
    std::size_t nhit_d;
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
    std::vector<std::size_t> *nhit_n = nullptr;
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
    std::vector<std::size_t> *nhit_mult = nullptr;
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
    std::size_t nhit_p;
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
    std::size_t nhit_d;
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
    std::vector<std::size_t> nhit_n;
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
    std::vector<std::size_t> nhit_mult;
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
        tin->SetBranchAddress("nhit_p", &varin.nhit_p);
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
        tin->SetBranchAddress("nhit_d", &varin.nhit_d);
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
        tin->SetBranchAddress("nhit_n", &varin.nhit_n);
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
        tin->SetBranchAddress("nhit_mult", &varin.nhit_mult);
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
        tout->SetBranchAddress("posx_p", &varout.posx_p);
        tout->SetBranchAddress("posy_p", &varout.posy_p);
        tout->SetBranchAddress("posz_p", &varout.posz_p);
        tout->SetBranchAddress("e_p", &varout.e_p);
        tout->SetBranchAddress("sec_p", &varout.sec_p);
        tout->SetBranchAddress("nsec_p", &varout.nsec_p);
        tout->SetBranchAddress("totq_p", &varout.totq_p);
        tout->SetBranchAddress("meanq_p", &varout.meanq_p);
        tout->SetBranchAddress("stdq_p", &varout.stdq_p);
        tout->SetBranchAddress("minq_p", &varout.minq_p);
        tout->SetBranchAddress("maxq_p", &varout.maxq_p);
        tout->SetBranchAddress("nhit_p", &varout.nhit_p);
        tout->SetBranchAddress("meant_p", &varout.meant_p);
        tout->SetBranchAddress("stdt_p", &varout.stdt_p);
        tout->SetBranchAddress("posx_d", &varout.posx_d);
        tout->SetBranchAddress("posy_d", &varout.posy_d);
        tout->SetBranchAddress("posz_d", &varout.posz_d);
        tout->SetBranchAddress("e_d", &varout.e_d);
        tout->SetBranchAddress("sec_d", &varout.sec_d);
        tout->SetBranchAddress("nsec_d", &varout.nsec_d);
        tout->SetBranchAddress("totq_d", &varout.totq_d);
        tout->SetBranchAddress("meanq_d", &varout.meanq_d);
        tout->SetBranchAddress("stdq_d", &varout.stdq_d);
        tout->SetBranchAddress("minq_d", &varout.minq_d);
        tout->SetBranchAddress("maxq_d", &varout.maxq_d);
        tout->SetBranchAddress("nhit_d", &varout.nhit_d);
        tout->SetBranchAddress("meant_d", &varout.meant_d);
        tout->SetBranchAddress("stdt_d", &varout.stdt_d);
        tout->SetBranchAddress("posx_n", &varout.posx_n);
        tout->SetBranchAddress("posy_n", &varout.posy_n);
        tout->SetBranchAddress("posz_n", &varout.posz_n);
        tout->SetBranchAddress("e_n", &varout.e_n);
        tout->SetBranchAddress("sec_n", &varout.sec_n);
        tout->SetBranchAddress("nsec_n", &varout.nsec_n);
        tout->SetBranchAddress("totq_n", &varout.totq_n);
        tout->SetBranchAddress("meanq_n", &varout.meanq_n);
        tout->SetBranchAddress("stdq_n", &varout.stdq_n);
        tout->SetBranchAddress("minq_n", &varout.minq_n);
        tout->SetBranchAddress("maxq_n", &varout.maxq_n);
        tout->SetBranchAddress("nhit_n", &varout.nhit_n);
        tout->SetBranchAddress("meant_n", &varout.meant_n);
        tout->SetBranchAddress("stdt_n", &varout.stdt_n);
        tout->SetBranchAddress("posx_mult", &varout.posx_mult);
        tout->SetBranchAddress("posy_mult", &varout.posy_mult);
        tout->SetBranchAddress("posz_mult", &varout.posz_mult);
        tout->SetBranchAddress("e_mult", &varout.e_mult);
        tout->SetBranchAddress("sec_mult", &varout.sec_mult);
        tout->SetBranchAddress("nsec_mult", &varout.nsec_mult);
        tout->SetBranchAddress("mult_type", &varout.mult_type);
        tout->SetBranchAddress("totq_mult", &varout.totq_mult);
        tout->SetBranchAddress("meanq_mult", &varout.meanq_mult);
        tout->SetBranchAddress("stdq_mult", &varout.stdq_mult);
        tout->SetBranchAddress("minq_mult", &varout.minq_mult);
        tout->SetBranchAddress("maxq_mult", &varout.maxq_mult);
        tout->SetBranchAddress("nhit_mult", &varout.nhit_mult);
        tout->SetBranchAddress("meant_mult", &varout.meant_mult);
        tout->SetBranchAddress("stdt_mult", &varout.stdt_mult);
        tout->SetBranchAddress("method_mu", &varout.method_mu);
        tout->SetBranchAddress("loc_mu", &varout.loc_mu);
        tout->SetBranchAddress("posx_mu", &varout.posx_mu);
        tout->SetBranchAddress("posy_mu", &varout.posy_mu);
        tout->SetBranchAddress("posz_mu", &varout.posz_mu);
        tout->SetBranchAddress("dirx_mu", &varout.dirx_mu);
        tout->SetBranchAddress("diry_mu", &varout.diry_mu);
        tout->SetBranchAddress("dirz_mu", &varout.dirz_mu);
        tout->SetBranchAddress("totq_mu", &varout.totq_mu);
        tout->SetBranchAddress("sec_mu", &varout.sec_mu);
        tout->SetBranchAddress("nsec_mu", &varout.nsec_mu);
        tout->SetBranchAddress("quality_mu", &varout.quality_mu);

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
            varout.nhit_p = varin.nhit_p;
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
            varout.nhit_d = varin.nhit_d;
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
            varout.nhit_n.assign(varin.nhit_n->end() - varin.e_n->size(), varin.nhit_n->end());
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
            varout.nhit_mult.assign(varin.nhit_mult->end() - varin.e_mult->size(), varin.nhit_mult->end());
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