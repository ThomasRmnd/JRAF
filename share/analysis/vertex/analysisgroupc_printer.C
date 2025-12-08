#include <chrono>
#include <iostream>
#include <set>
#include <string>

#include <TFile.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TVector3.h>

TTimeStamp operator-(const TTimeStamp& lhs, const TTimeStamp& rhs) {
    return TTimeStamp{
        lhs.GetSec() - rhs.GetSec(), 
        lhs.GetNanoSec() - rhs.GetNanoSec()
    };
}

TTimeStamp operator+(const TTimeStamp& lhs, const TTimeStamp& rhs) {
    return TTimeStamp{
        lhs.GetSec() + rhs.GetSec(), 
        lhs.GetNanoSec() + rhs.GetNanoSec()
    };
}

class AnalysisBase {

public:

    std::string name;

    int run_id;

    double posx_p;
    double posy_p;
    double posz_p;
    double e_p;
    double totq_p;
    time_t sec_p;
    int nsec_p;

    double posx_d;
    double posy_d;
    double posz_d;
    double e_d;
    double totq_d;
    time_t sec_d;
    int nsec_d;

    AnalysisBase(const std::string& name_) : name{name_} {}

    virtual ~AnalysisBase() = default;

    virtual bool bind(TFile* file, TTree** tree) {
        *tree = file->Get<TTree>(name.c_str());
        if (!*tree) {
            std::cerr << "Tree " << name << " not found in file\n";
            return false;
        }
        // (*tree)->SetBranchAddress("run_id", &run_id);

        (*tree)->SetBranchAddress("posx_p", &posx_p);
        (*tree)->SetBranchAddress("posy_p", &posy_p);
        (*tree)->SetBranchAddress("posz_p", &posz_p);
        (*tree)->SetBranchAddress("e_p", &e_p);
        (*tree)->SetBranchAddress("totq_p", &totq_p);
        (*tree)->SetBranchAddress("sec_p", &sec_p);
        (*tree)->SetBranchAddress("nsec_p", &nsec_p);

        (*tree)->SetBranchAddress("posx_d", &posx_d);
        (*tree)->SetBranchAddress("posy_d", &posy_d);
        (*tree)->SetBranchAddress("posz_d", &posz_d);
        (*tree)->SetBranchAddress("e_d", &e_d);
        (*tree)->SetBranchAddress("totq_d", &totq_d);
        (*tree)->SetBranchAddress("sec_d", &sec_d);
        (*tree)->SetBranchAddress("nsec_d", &nsec_d);

        return true;
    }

    virtual bool selection() = 0;

    virtual void print() {
        std::cout << "Prompt: " /* (" << posx_p << ", " << posy_p << ", " << posz_p << "), */ "E = " << e_p << ", Q = " << totq_p << ", Time = " << TTimeStamp(sec_p, nsec_p) << '\n';
        std::cout << "Delayed: " /* (" << posx_d << ", " << posy_d << ", " << posz_d << "), */ "E = " << e_d << ", Q = " << totq_d << ", Time = " << TTimeStamp(sec_d, nsec_d) << '\n';
    }

};

class IBDAnalysis : public AnalysisBase {

public:

    IBDAnalysis(const std::string& suffix) : AnalysisBase{"IBDAnalysis" + suffix} {}

    ~IBDAnalysis() override = default;

    std::vector<double>* posx_n = nullptr;
    std::vector<double>* posy_n = nullptr;
    std::vector<double>* posz_n = nullptr;
    std::vector<double>* e_n = nullptr;
    std::vector<double>* totq_n = nullptr;
    std::vector<time_t>* sec_n = nullptr;
    std::vector<int>* nsec_n = nullptr;

    std::vector<double>* posx_mult = nullptr;
    std::vector<double>* posy_mult = nullptr;
    std::vector<double>* posz_mult = nullptr;
    std::vector<double>* e_mult = nullptr;
    std::vector<double>* totq_mult = nullptr;
    std::vector<time_t>* sec_mult = nullptr;
    std::vector<int>* nsec_mult = nullptr;
    std::vector<int>* mult_type = nullptr;

    std::vector<std::string>* method_mu = nullptr;
    std::vector<int>* loc_mu = nullptr;
    std::vector<double>* posx_mu = nullptr;
    std::vector<double>* posy_mu = nullptr;
    std::vector<double>* posz_mu = nullptr;
    std::vector<double>* dirx_mu = nullptr;
    std::vector<double>* diry_mu = nullptr;
    std::vector<double>* dirz_mu = nullptr;
    std::vector<double>* totq_mu = nullptr;
    std::vector<time_t>* sec_mu = nullptr;
    std::vector<int>* nsec_mu = nullptr;
    std::vector<double>* quality_mu = nullptr;

    bool bind(TFile* file, TTree** tree) override {
        if (!AnalysisBase::bind(file, tree)) return false;
        
        (*tree)->SetBranchAddress("posx_n", &posx_n);
        (*tree)->SetBranchAddress("posy_n", &posy_n);
        (*tree)->SetBranchAddress("posz_n", &posz_n);
        (*tree)->SetBranchAddress("e_n", &e_n);
        (*tree)->SetBranchAddress("totq_n", &totq_n);
        (*tree)->SetBranchAddress("sec_n", &sec_n);
        (*tree)->SetBranchAddress("nsec_n", &nsec_n);

        (*tree)->SetBranchAddress("posx_mult", &posx_mult);
        (*tree)->SetBranchAddress("posy_mult", &posy_mult);
        (*tree)->SetBranchAddress("posz_mult", &posz_mult);
        (*tree)->SetBranchAddress("e_mult", &e_mult);
        (*tree)->SetBranchAddress("totq_mult", &totq_mult);
        (*tree)->SetBranchAddress("sec_mult", &sec_mult);
        (*tree)->SetBranchAddress("nsec_mult", &nsec_mult);
        (*tree)->SetBranchAddress("mult_type", &mult_type);

        (*tree)->SetBranchAddress("method_mu", &method_mu);
        (*tree)->SetBranchAddress("loc_mu", &loc_mu);
        (*tree)->SetBranchAddress("posx_mu", &posx_mu);
        (*tree)->SetBranchAddress("posy_mu", &posy_mu);
        (*tree)->SetBranchAddress("posz_mu", &posz_mu);
        (*tree)->SetBranchAddress("dirx_mu", &dirx_mu);
        (*tree)->SetBranchAddress("diry_mu", &diry_mu);
        (*tree)->SetBranchAddress("dirz_mu", &dirz_mu);
        (*tree)->SetBranchAddress("totq_mu", &totq_mu);
        (*tree)->SetBranchAddress("sec_mu", &sec_mu);
        (*tree)->SetBranchAddress("nsec_mu", &nsec_mu);
        (*tree)->SetBranchAddress("quality_mu", &quality_mu);

        return true;
    }

    bool selection() override {
        TTimeStamp ts_p{sec_p, nsec_p};
        TTimeStamp ts_d{sec_d, nsec_d};
        TVector3 pos_p{posx_p, posy_p, posz_p};
        TVector3 pos_d{posx_d, posy_d, posz_d};
        if (pos_p.Mag() > 16500.0) return false;
        if ((posz_p < -15500.0 || 15500 < posz_p) && std::sqrt(posx_p * posx_p + posy_p * posy_p) < 3000.0) return false;
        if (e_p < 0.7 || 12.0 < e_p) return false;
        if (e_d < 2.0 || 2.5 < e_d) return false;

        std::size_t nb_neutron_veto = 0ul;
        for (std::size_t k = 0ul; k < e_n->size(); ++k) {
            if (e_n->operator[](k) < 1.5 || 20.0 < e_n->operator[](k)) continue;
            TTimeStamp ts_n{sec_n->operator[](k), nsec_n->operator[](k)};
            TVector3 pos_n{posx_n->operator[](k), posy_n->operator[](k), posz_n->operator[](k)};
            if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_p - ts_n < TTimeStamp{0, 20000} || TTimeStamp{0, 1200000000} < ts_p - ts_n) continue;
            if (ts_d - ts_n < TTimeStamp{0, 20000} || TTimeStamp{0, 1200000000} < ts_d - ts_n) continue;
            ++nb_neutron_veto;
        }

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < e_mult->size(); ++k) {
            if (e_mult->operator[](k) < 2.0 || 12.0 < e_mult->operator[](k)) continue;
            TTimeStamp ts_mult{sec_mult->operator[](k), nsec_mult->operator[](k)};
            TVector3 pos_mult{posx_mult->operator[](k), posy_mult->operator[](k), posz_mult->operator[](k)};
            // if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_p - ts_mult < TTimeStamp{0, 1000000} || TTimeStamp{0, 0} < ts_p - ts_mult) continue;
            if (ts_d - ts_mult < TTimeStamp{0, 0} || TTimeStamp{0, 1000000} < ts_d - ts_mult) continue;
            ++nb_multu_veto;
        }

        return (nb_neutron_veto == 0ul && nb_multu_veto == 0ul);
    }

    void print() override {
        TTimeStamp ts_p{sec_p, nsec_p};
        TTimeStamp ts_d{sec_d, nsec_d};
        TVector3 pos_p{posx_p, posy_p, posz_p};
        TVector3 pos_d{posx_d, posy_d, posz_d};

        std::string method = "";
        double d_mu2p = std::numeric_limits<double>::infinity();
        double t_mu2p = std::numeric_limits<double>::infinity();
        double d_mu2d = std::numeric_limits<double>::infinity();
        double t_mu2d = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0ul; k < method_mu->size(); ++k) {
            TTimeStamp ts_mu{sec_mu->operator[](k), nsec_mu->operator[](k)};
            TVector3 pos_mu{posx_mu->operator[](k), posy_mu->operator[](k), posz_mu->operator[](k)};
            TVector3 dir_mu{dirx_mu->operator[](k), diry_mu->operator[](k), dirz_mu->operator[](k)};
            double tmp_d_mu2p = dir_mu.Cross(pos_p - pos_mu).Mag();
            double tmp_t_mu2p = static_cast<double>(ts_p - ts_mu);
            double tmp_d_mu2d = dir_mu.Cross(pos_d - pos_mu).Mag();
            double tmp_t_mu2d = static_cast<double>(ts_d - ts_mu);
            if (tmp_d_mu2p < d_mu2p && 0.0 < tmp_t_mu2p && tmp_t_mu2p < 1.2) {
                method = method_mu->operator[](k);
                d_mu2p = tmp_d_mu2p;
                t_mu2p = tmp_t_mu2p;
                d_mu2d = tmp_d_mu2d;
                t_mu2d = tmp_t_mu2d;
            }
        }

        if (ts_p < TTimeStamp(2025, 8, 30, 0, 0, 0) || TTimeStamp(2025, 8, 31, 0, 0, 0) < ts_p) return; 

        // std::cout << "Prompt: (" << posx_p << ", " << posy_p << ", " << posz_p << "), E = " << e_p << ", Q = " << totq_p << ", Time = " << TTimeStamp(sec_p, nsec_p) << '\n';
        // std::cout << "Delayed: (" << posx_d << ", " << posy_d << ", " << posz_d << "), E = " << e_d << ", Q = " << totq_d << ", Time = " << TTimeStamp(sec_d, nsec_d) << '\n';
        // std::cout << "Number of Neutron Veto associated: " << nb_neutron_veto << '\n';
        // std::cout << "Number of Multiplicity Veto associated: " << nb_multu_veto << '\n';
        // std::cout << "Muon (" << method << "): d_mu2p = " << d_mu2p << ", t_mu2p = " << t_mu2p << ", d_mu2d = " << d_mu2d << ", t_mu2d = " << t_mu2d << '\n';

        if (selection()) {
            std::cout << "Prompt: E = " << e_p << ", Q = " << totq_p << ", Time = " << TTimeStamp{sec_p, nsec_p} << '\n';
            std::cout << "Delayed: E = " << e_d << ", Q = " << totq_d << ", Time = " << TTimeStamp{sec_d, nsec_d} << '\n';
        }
    }

};

struct ThomasIBD {

    int run_id;

    TTimeStamp ts_p;
    double e_p;
    double totq_p;

    TTimeStamp ts_d;
    double e_d;
    double totq_d;

};

bool operator<(const ThomasIBD& lhs, const ThomasIBD& rhs) {
    return lhs.ts_p < rhs.ts_p;
}

void print_all_entries(TFile* file, TTree** tree, AnalysisBase* analysis) {
    if (!analysis->bind(file, tree)) return;
    std::set<ThomasIBD> ibds;
    std::cout << "=== Analysis: " << analysis->name << " (Total Entries: " << (*tree)->GetEntries() << ") ===\n";
    for (long k = 0; k < (*tree)->GetEntries(); ++k) {
        (*tree)->GetEntry(k);
        // analysis->print();
        if (!analysis->selection()) continue;
        ThomasIBD ibd;
        ibd.run_id = 0;
        ibd.ts_p = TTimeStamp{analysis->sec_p, analysis->nsec_p};
        ibd.e_p = analysis->e_p;
        ibd.totq_p = analysis->totq_p;
        ibd.ts_d = TTimeStamp{analysis->sec_d, analysis->nsec_d};
        ibd.e_d = analysis->e_d;
        ibd.totq_d = analysis->totq_d;
        ibds.insert(ibd);
    }

    for (std::set<ThomasIBD>::const_iterator it = ibds.begin(); it != ibds.end(); ++it) {
        std::cout << "Prompt: " << "E = " << it->e_p << ", Q = " << it->totq_p << ", Time = " << it->ts_p << '\n';
        std::cout << "Delayed: " << "E = " << it->e_d << ", Q = " << it->totq_d << ", Time = " << it->ts_d << '\n';
    }
}


struct VanessaIBD {

    int run_id;

    TTimeStamp ts_p;
    double e_p;
    double totq_p;

    TTimeStamp ts_d;
    double e_d;
    double totq_d;

};

bool operator<(const VanessaIBD& lhs, const VanessaIBD& rhs) {
    return lhs.ts_p < rhs.ts_p;
}

void analyze_vanessa_result(TFile* file, TTree* tree) {
    tree = file->Get<TTree>("events");
    std::set<VanessaIBD> ibds;
    double run_id;
    double t_p, t_d;
    double e_p, e_d;
    double totq_p, totq_d;
    tree->SetBranchAddress("run_number", &run_id);
    tree->SetBranchAddress("time_p_ns", &t_p);
    tree->SetBranchAddress("time_d_ns", &t_d);
    tree->SetBranchAddress("energy_p_omilrec", &e_p);
    tree->SetBranchAddress("energy_d_omilrec", &e_d);
    tree->SetBranchAddress("NPE_p", &totq_p);
    tree->SetBranchAddress("NPE_d", &totq_d);
    for (long k = 0l; k < tree->GetEntries(); ++k) {
        tree->GetEntry(k);
        time_t sec_p = static_cast<time_t>(t_p / 1.0e9);
        int nsec_p = static_cast<int>(t_p - static_cast<double>(sec_p) * 1.0e9);
        time_t sec_d = static_cast<time_t>(t_d / 1.0e9);
        int nsec_d = static_cast<int>(t_d - static_cast<double>(sec_d) * 1.0e9);
        VanessaIBD ibd;
        ibd.run_id = run_id;
        ibd.ts_p = TTimeStamp{sec_p, nsec_p};
        ibd.e_p = e_p;
        ibd.totq_p = totq_p;
        ibd.ts_d = TTimeStamp{sec_d, nsec_d};
        ibd.e_d = e_d;
        ibd.totq_d = totq_d;
        ibds.insert(ibd);
    }
    for (std::set<VanessaIBD>::const_iterator it = ibds.begin(); it != ibds.end(); ++it) {
        std::cout << "(RUN = " << it->run_id << ") Prompt: " << "E = " << it->e_p << ", Q = " << it->totq_p << ", Time = " << it->ts_p << '\n';
        std::cout << "(RUN = " << it->run_id << ") Delayed: " << "E = " << it->e_d << ", Q = " << it->totq_d << ", Time = " << it->ts_d << '\n';
    }
}


void analysisgroupc_printer(const std::string& filename, const std::string& suffix) {
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open file: " << filename << '\n';
        return;
    }
    TTree* tree = nullptr;

    // analyze_vanessa_result(file, tree);

    /*

    IBDWithCylindricalCut ibd_with_cylindrical_cut_3m{"_3m" + suffix};
    print_all_entries(file, &tree, &ibd_with_cylindrical_cut_3m);

    IBDWithCylindricalCut ibd_with_cylindrical_cut_5m{"_5m" + suffix};
    print_all_entries(file, &tree, &ibd_with_cylindrical_cut_5m);

    */

    IBDAnalysis ibd_analysis(suffix);
    print_all_entries(file, &tree, &ibd_analysis);

    /*

    MultiplicityWindowCut multiplicity_window_cut(suffix);
    print_all_entries(file, &tree, &multiplicity_window_cut);

    CosmoRateWithNeutronAnalysis cosmo_rate_with_neutron_study(suffix);
    print_all_entries(file, &tree, &cosmo_rate_with_neutron_study);

    TtCosmoStudy tt_cosmo_study_sig{"_sig" + suffix};
    print_all_entries(file, &tree, &tt_cosmo_study_sig);

    TtCosmoStudy tt_cosmo_study_bkg{"_bkg" + suffix};
    print_all_entries(file, &tree, &tt_cosmo_study_bkg);

    CdWpCosmoStudy cdwp_cosmo_study_sig{"_sig" + suffix};
    print_all_entries(file, &tree, &cdwp_cosmo_study_sig);

    CdWpCosmoStudy cdwp_cosmo_study_bkg{"_bkg" + suffix};
    print_all_entries(file, &tree, &cdwp_cosmo_study_bkg);

    CdWpCosmoStudy cdwp_all_cosmo_study{"_All" + suffix};
    print_all_entries(file, &tree, &cdwp_all_cosmo_study);

    */

    file->Close();
    delete file;
}