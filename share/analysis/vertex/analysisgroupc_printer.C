#include <algorithm>
#include <chrono>
#include <iostream>
#include <set>
#include <string>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
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

struct Vertex {

    TVector3 pos;
    TTimeStamp ts;
    double e;
    double q;

};

inline bool operator<(const Vertex& lhs, const Vertex& rhs) {
    return lhs.ts < rhs.ts;
}

struct IBD {

    Vertex prompt;
    Vertex delayed;

};

inline bool operator<(const IBD& lhs, const IBD& rhs) {
    return lhs.prompt.ts < rhs.prompt.ts;
}

class AnalysisBase {

public:

    std::string name;

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

    AnalysisBase(const std::string& name_) : name{name_} {}

    virtual ~AnalysisBase() = default;

    virtual TChain* retrieve(const std::string& filename) {
        TChain* chain = new TChain(name.c_str());
        if (!chain) {
            std::cerr << "Cannot create TChain " << name << '\n';
            return nullptr;
        }
        chain->Add(filename.c_str());

        chain->SetBranchAddress("posx_p", &posx_p);
        chain->SetBranchAddress("posy_p", &posy_p);
        chain->SetBranchAddress("posz_p", &posz_p);
        chain->SetBranchAddress("e_p", &e_p);
        chain->SetBranchAddress("sec_p", &sec_p);
        chain->SetBranchAddress("nsec_p", &nsec_p);

        chain->SetBranchAddress("totq_p", &totq_p);
        chain->SetBranchAddress("meanq_p", &meanq_p);
        chain->SetBranchAddress("stdq_p", &stdq_p);
        chain->SetBranchAddress("minq_p", &minq_p);
        chain->SetBranchAddress("maxq_p", &maxq_p);
        chain->SetBranchAddress("nhit_p", &nhit_p);
        chain->SetBranchAddress("meant_p", &meant_p);
        chain->SetBranchAddress("stdt_p", &stdt_p);

        chain->SetBranchAddress("posx_d", &posx_d);
        chain->SetBranchAddress("posy_d", &posy_d);
        chain->SetBranchAddress("posz_d", &posz_d);
        chain->SetBranchAddress("e_d", &e_d);
        chain->SetBranchAddress("sec_d", &sec_d);
        chain->SetBranchAddress("nsec_d", &nsec_d);

        chain->SetBranchAddress("totq_d", &totq_d);
        chain->SetBranchAddress("meanq_d", &meanq_d);
        chain->SetBranchAddress("stdq_d", &stdq_d);
        chain->SetBranchAddress("minq_d", &minq_d);
        chain->SetBranchAddress("maxq_d", &maxq_d);
        chain->SetBranchAddress("nhit_d", &nhit_d);
        chain->SetBranchAddress("meant_d", &meant_d);
        chain->SetBranchAddress("stdt_d", &stdt_d);

        return chain;
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

    TChain* retrieve(const std::string& filename) override {
        TChain* chain = AnalysisBase::retrieve(filename);
        if (!chain) return nullptr;
        
        chain->SetBranchAddress("posx_n", &posx_n);
        chain->SetBranchAddress("posy_n", &posy_n);
        chain->SetBranchAddress("posz_n", &posz_n);
        chain->SetBranchAddress("e_n", &e_n);
        chain->SetBranchAddress("totq_n", &totq_n);
        chain->SetBranchAddress("sec_n", &sec_n);
        chain->SetBranchAddress("nsec_n", &nsec_n);

        chain->SetBranchAddress("posx_mult", &posx_mult);
        chain->SetBranchAddress("posy_mult", &posy_mult);
        chain->SetBranchAddress("posz_mult", &posz_mult);
        chain->SetBranchAddress("e_mult", &e_mult);
        chain->SetBranchAddress("totq_mult", &totq_mult);
        chain->SetBranchAddress("sec_mult", &sec_mult);
        chain->SetBranchAddress("nsec_mult", &nsec_mult);
        chain->SetBranchAddress("mult_type", &mult_type);

        chain->SetBranchAddress("method_mu", &method_mu);
        chain->SetBranchAddress("loc_mu", &loc_mu);
        chain->SetBranchAddress("posx_mu", &posx_mu);
        chain->SetBranchAddress("posy_mu", &posy_mu);
        chain->SetBranchAddress("posz_mu", &posz_mu);
        chain->SetBranchAddress("dirx_mu", &dirx_mu);
        chain->SetBranchAddress("diry_mu", &diry_mu);
        chain->SetBranchAddress("dirz_mu", &dirz_mu);
        chain->SetBranchAddress("totq_mu", &totq_mu);
        chain->SetBranchAddress("sec_mu", &sec_mu);
        chain->SetBranchAddress("nsec_mu", &nsec_mu);
        chain->SetBranchAddress("quality_mu", &quality_mu);

        return chain;
    }

    bool selection() override {
        if (stdt_p > 200.0 || stdt_d > 200.0) return false;
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
            if (ts_mult < ts_p - TTimeStamp{0, 1000000} || ts_d + TTimeStamp{0, 1000000} < ts_mult) continue;
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

std::vector<IBD> get_all_ibd(const std::string& filename, AnalysisBase* analysis) {
    TChain* chain = analysis->retrieve(filename);
    if (!chain) return {};
    std::set<IBD> ibds_ordered;
    std::cout << "=== Analysis: " << analysis->name << " (Total Entries: " << chain->GetEntries() << ") ===\n";
    for (long k = 0; k < chain->GetEntries(); ++k) {
        chain->GetEntry(k);
        // analysis->print();
        if (!analysis->selection()) continue;
        IBD ibd;
        ibd.prompt.pos = TVector3{analysis->posx_p, analysis->posy_p, analysis->posz_p};
        ibd.prompt.ts = TTimeStamp{analysis->sec_p, analysis->nsec_p};
        ibd.prompt.e = analysis->e_p;
        ibd.prompt.q = analysis->totq_p;
        ibd.delayed.pos = TVector3{analysis->posx_d, analysis->posy_d, analysis->posz_d};
        ibd.delayed.ts = TTimeStamp{analysis->sec_d, analysis->nsec_d};
        ibd.delayed.e = analysis->e_d;
        ibd.delayed.q = analysis->totq_d;
        ibds_ordered.insert(ibd);
    }

    std::vector<IBD> ibds;
    ibds.reserve(ibds_ordered.size());
    for (std::set<IBD>::const_iterator it = ibds_ordered.begin(); it != ibds_ordered.end(); ++it) {
        ibds.push_back(*it);
    }
    return ibds;
}

void print_all_entries(const std::string& filename, AnalysisBase* analysis) {
    TChain* chain = analysis->retrieve(filename);
    if (!chain) return;
    std::set<IBD> ibds;
    std::cout << "=== Analysis: " << analysis->name << " (Total Entries: " << chain->GetEntries() << ") ===\n";
    for (long k = 0; k < chain->GetEntries(); ++k) {
        chain->GetEntry(k);
        // analysis->print();
        if (!analysis->selection()) continue;
        IBD ibd;
        ibd.prompt.pos = TVector3{analysis->posx_p, analysis->posy_p, analysis->posz_p};
        ibd.prompt.ts = TTimeStamp{analysis->sec_p, analysis->nsec_p};
        ibd.prompt.e = analysis->e_p;
        ibd.prompt.q = analysis->totq_p;
        ibd.delayed.pos = TVector3{analysis->posx_d, analysis->posy_d, analysis->posz_d};
        ibd.delayed.ts = TTimeStamp{analysis->sec_d, analysis->nsec_d};
        ibd.delayed.e = analysis->e_d;
        ibd.delayed.q = analysis->totq_d;
        ibds.insert(ibd);
    }

    for (std::set<IBD>::const_iterator it = ibds.begin(); it != ibds.end(); ++it) {
        std::cout << "Prompt: " << "E = " << it->prompt.e << ", Q = " << it->prompt.q << ", Time = " << it->prompt.ts << '\n';
        std::cout << "Delayed: " << "E = " << it->delayed.e << ", Q = " << it->delayed.q << ", Time = " << it->delayed.ts << '\n';
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

void compare_with_vanessa(const std::string& filename, IBDAnalysis* analysis) {
    TFile* vanessa_file = TFile::Open("/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/ibd/summary/ReProd25C/IBD_all_reprodC.root", "READ");
    TTree* vanessa_tree = vanessa_file->Get<TTree>("events");
    std::set<VanessaIBD> vanessa_ibds;
    double run_id;
    double t_p, t_d;
    double e_p, e_d;
    double totq_p, totq_d;
    vanessa_tree->SetBranchAddress("run_number", &run_id);
    vanessa_tree->SetBranchAddress("time_p_ns", &t_p);
    vanessa_tree->SetBranchAddress("time_d_ns", &t_d);
    vanessa_tree->SetBranchAddress("energy_p_omilrec", &e_p);
    vanessa_tree->SetBranchAddress("energy_d_omilrec", &e_d);
    vanessa_tree->SetBranchAddress("NPE_p", &totq_p);
    vanessa_tree->SetBranchAddress("NPE_d", &totq_d);
    for (long k = 0l; k < vanessa_tree->GetEntries(); ++k) {
        vanessa_tree->GetEntry(k);
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
        vanessa_ibds.insert(ibd);
    }

    TChain* chain = analysis->retrieve(filename);
    if (!chain) return;
    std::set<IBD> ibds;
    std::cout << "=== Analysis: " << analysis->name << " (Total Entries: " << chain->GetEntries() << ") ===\n";
    for (long k = 0; k < chain->GetEntries(); ++k) {
        chain->GetEntry(k);
        std::set<VanessaIBD>::const_iterator it = std::find_if(
            vanessa_ibds.begin(),
            vanessa_ibds.end(),
            [&](const VanessaIBD& vanessa_ibd) {
                return (
                    vanessa_ibd.ts_p.GetSec() == analysis->sec_p &&
                    vanessa_ibd.ts_d.GetSec() == analysis->sec_d
                );
            }
        );
        bool is_only_in_vanessa = (it != vanessa_ibds.end() && !analysis->selection());
        bool is_only_in_analysis = (it == vanessa_ibds.end() && analysis->selection());
        if (!is_only_in_vanessa && !is_only_in_analysis) continue;

        TTimeStamp ts_p{analysis->sec_p, analysis->nsec_p};
        TTimeStamp ts_d{analysis->sec_d, analysis->nsec_d};
        TVector3 pos_p{analysis->posx_p, analysis->posy_p, analysis->posz_p};
        TVector3 pos_d{analysis->posx_d, analysis->posy_d, analysis->posz_d};

        std::size_t nb_neutron_veto = 0ul;
        for (std::size_t k = 0ul; k < analysis->e_n->size(); ++k) {
            if (analysis->e_n->operator[](k) < 1.5 || 20.0 < analysis->e_n->operator[](k)) continue;
            TTimeStamp ts_n{analysis->sec_n->operator[](k), analysis->nsec_n->operator[](k)};
            TVector3 pos_n{analysis->posx_n->operator[](k), analysis->posy_n->operator[](k), analysis->posz_n->operator[](k)};
            if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_p - ts_n < TTimeStamp{0, 20000} || TTimeStamp{0, 1200000000} < ts_p - ts_n) continue;
            if (ts_d - ts_n < TTimeStamp{0, 20000} || TTimeStamp{0, 1200000000} < ts_d - ts_n) continue;
            ++nb_neutron_veto;
        }

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < analysis->e_mult->size(); ++k) {
            if (analysis->e_mult->operator[](k) < 2.0 || 12.0 < analysis->e_mult->operator[](k)) continue;
            TTimeStamp ts_mult{analysis->sec_mult->operator[](k), analysis->nsec_mult->operator[](k)};
            TVector3 pos_mult{analysis->posx_mult->operator[](k), analysis->posy_mult->operator[](k), analysis->posz_mult->operator[](k)};
            // if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_mult < ts_p - TTimeStamp{0, 1000000} || ts_d + TTimeStamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }

        std::string method = "";
        double d_mu2p = std::numeric_limits<double>::infinity();
        double t_mu2p = std::numeric_limits<double>::infinity();
        double d_mu2d = std::numeric_limits<double>::infinity();
        double t_mu2d = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0ul; k < analysis->method_mu->size(); ++k) {
            TTimeStamp ts_mu{analysis->sec_mu->operator[](k), analysis->nsec_mu->operator[](k)};
            TVector3 pos_mu{analysis->posx_mu->operator[](k), analysis->posy_mu->operator[](k), analysis->posz_mu->operator[](k)};
            TVector3 dir_mu{analysis->dirx_mu->operator[](k), analysis->diry_mu->operator[](k), analysis->dirz_mu->operator[](k)};
            double tmp_d_mu2p = dir_mu.Cross(pos_p - pos_mu).Mag();
            double tmp_t_mu2p = static_cast<double>(ts_p - ts_mu);
            double tmp_d_mu2d = dir_mu.Cross(pos_d - pos_mu).Mag();
            double tmp_t_mu2d = static_cast<double>(ts_d - ts_mu);
            if (tmp_d_mu2p < d_mu2p && 0.0 < tmp_t_mu2p && tmp_t_mu2p < 1.2) {
                method = analysis->method_mu->operator[](k);
                d_mu2p = tmp_d_mu2p;
                t_mu2p = tmp_t_mu2p;
                d_mu2d = tmp_d_mu2d;
                t_mu2d = tmp_t_mu2d;
            }
        }
        
        std::cout << (is_only_in_vanessa ? "[Only in Vanessa] ======================================" : "");
        std::cout << (is_only_in_analysis ? "[Only in Analysis] ======================================" : "");
        std::cout << '\n';
        std::cout << "Prompt: (" << analysis->posx_p << ", " << analysis->posy_p << ", " << analysis->posz_p << "), E = " << analysis->e_p << ", Q = " << analysis->totq_p << ", Time = " << TTimeStamp(analysis->sec_p, analysis->nsec_p) << '\n';
        std::cout << "Delayed: (" << analysis->posx_d << ", " << analysis->posy_d << ", " << analysis->posz_d << "), E = " << analysis->e_d << ", Q = " << analysis->totq_d << ", Time = " << TTimeStamp(analysis->sec_d, analysis->nsec_d) << '\n';
        std::cout << "Number of Neutron Veto associated: " << nb_neutron_veto << '\n';
        std::cout << "Number of Multiplicity Veto associated: " << nb_multu_veto << '\n';

        // print neutrons
        for (std::size_t k = 0ul; k < analysis->e_n->size(); ++k) {
            TTimeStamp ts_n{analysis->sec_n->operator[](k), analysis->nsec_n->operator[](k)};
            TVector3 pos_n{analysis->posx_n->operator[](k), analysis->posy_n->operator[](k), analysis->posz_n->operator[](k)};
            double e_n = analysis->e_n->operator[](k);
            double totq_n = analysis->totq_n->operator[](k);
            std::cout << "  Neutron: (" << pos_n.X() << ", " << pos_n.Y() << ", " << pos_n.Z() << "), E = " << e_n << ", Q = " << totq_n << ", Time = " << ts_n << '\n';
        }

        // print mults
        for (std::size_t k = 0ul; k < analysis->e_mult->size(); ++k) {
            TTimeStamp ts_mult{analysis->sec_mult->operator[](k), analysis->nsec_mult->operator[](k)};
            TVector3 pos_mult{analysis->posx_mult->operator[](k), analysis->posy_mult->operator[](k), analysis->posz_mult->operator[](k)};
            double e_mult = analysis->e_mult->operator[](k);
            double totq_mult = analysis->totq_mult->operator[](k);
            std::cout << "  Mult: (" << pos_mult.X() << ", " << pos_mult.Y() << ", " << pos_mult.Z() << "), E = " << e_mult << ", Q = " << totq_mult << ", Time = " << ts_mult << '\n';
        }

        // print muons
        for (std::size_t k = 0ul; k < analysis->method_mu->size(); ++k) {
            TTimeStamp ts_mu{analysis->sec_mu->operator[](k), analysis->nsec_mu->operator[](k)};
            TVector3 pos_mu{analysis->posx_mu->operator[](k), analysis->posy_mu->operator[](k), analysis->posz_mu->operator[](k)};
            TVector3 dir_mu{analysis->dirx_mu->operator[](k), analysis->diry_mu->operator[](k), analysis->dirz_mu->operator[](k)};
            double totq_mu = analysis->totq_mu->operator[](k);
            double quality_mu = analysis->quality_mu->operator[](k);
            std::cout << "  Muon (" << analysis->method_mu->operator[](k) << "): (" << pos_mu.X() << ", " << pos_mu.Y() << ", " << pos_mu.Z() << "), Dir = (" << dir_mu.X() << ", " << dir_mu.Y() << ", " << dir_mu.Z() << "), Q = " << totq_mu << ", Quality = " << quality_mu << ", Time = " << ts_mu << '\n';
        }
    }
}

std::vector<double> linspace_cpp(double start, double stop, int num) {
    if (num <= 1) {
        if (num == 1) return {start};
        std::cerr << "Warning: linspace_cpp requires num >= 1. Returning empty vector\n";
        return {};
    }
    double step = (stop - start) / (num - 1);
    std::vector<double> result;
    for (int i = 0; i < num; ++i) {
        double value = start + static_cast<double>(i) * step;
        if (i == num - 1) {
            result.push_back(stop);
        } 
        else {
            result.push_back(value);
        }
    }
    
    return result;
}

#undef PRINT_ALL_ENTRIES
#undef COMPARE_WITH_VANESSA
#define GET_ALL_IBD

void analysisgroupc_printer(const std::string& filename, const std::string& suffix) {
    // analyze_vanessa_result(file, tree);

    IBDAnalysis ibd_analysis(suffix);
#ifdef PRINT_ALL_ENTRIES
    print_all_entries(filename, &ibd_analysis);
#endif 
#ifdef COMPARE_WITH_VANESSA
    compare_with_vanessa(filename, &ibd_analysis);
#endif
#ifdef GET_ALL_IBD
    std::vector<IBD> ibds = get_all_ibd(filename, &ibd_analysis);


    double e_p_min = 0.7;
    double e_p_max = 12.0;
    double e_p_width = 0.20;
    int e_p_nbin = std::round((e_p_max - e_p_min) / e_p_width) + 1;
    std::vector<double> e_p_bins = linspace_cpp(e_p_min, e_p_max, e_p_nbin);
    TH1D* h_e_p = new TH1D("h_e_p", "Prompt energy", e_p_bins.size() - 1, e_p_bins.data());

    double e_d_min = 2.0;
    double e_d_max = 2.5;
    double e_d_width = 0.02;
    int e_d_nbin = std::round((e_d_max - e_d_min) / e_d_width) + 1;
    std::vector<double> e_d_bins = linspace_cpp(e_d_min, e_d_max, e_d_nbin);
    TH1D* h_e_d = new TH1D("h_e_d", "Delayed energy", e_d_bins.size() - 1, e_d_bins.data());

    double e_dt_min = 0.0;
    double e_dt_max = 1.0;
    double e_dt_width = 0.025;
    int e_dt_nbin = std::round((e_dt_max - e_dt_min) / e_dt_width) + 1;
    std::vector<double> e_dt_bins = linspace_cpp(e_dt_min, e_dt_max, e_dt_nbin);
    TH1D* h_dt = new TH1D("h_dt", "Prompt-Delayed time difference", e_dt_bins.size() - 1, e_dt_bins.data());

    double e_dr_min = 0.0;
    double e_dr_max = 1.5;
    double e_dr_width = 0.05;
    int e_dr_nbin = std::round((e_dr_max - e_dr_min) / e_dr_width) + 1;
    std::vector<double> e_dr_bins = linspace_cpp(e_dr_min, e_dr_max, e_dr_nbin);
    TH1D* h_dr = new TH1D("h_dr", "Prompt-Delayed distance", e_dr_bins.size() - 1, e_dr_bins.data());

    double rho_min = 0.0;
    double rho_max = 17.7 * 17.7;
    int rho_nbin = 51;
    double z_min = -20.0;
    double z_max = 20.0;
    int z_nbin = 51;
    std::vector<double> rho_bins = linspace_cpp(rho_min, rho_max, rho_nbin);
    std::vector<double> z_bins = linspace_cpp(z_min, z_max, z_nbin);
    TH2D* h_rho_z_p = new TH2D("h_rho_z_p", "Prompt vertex distribution", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
    TH2D* h_rho_z_d = new TH2D("h_rho_z_d", "Delayed vertex distribution", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());

    for (const IBD& ibd : ibds) {
        h_e_p->Fill(ibd.prompt.e);
        h_e_d->Fill(ibd.delayed.e);
        h_dt->Fill((ibd.delayed.ts - ibd.prompt.ts) * 1000.0);
        h_dr->Fill((ibd.delayed.pos - ibd.prompt.pos).Mag() / 1000.0);
        h_rho_z_p->Fill((ibd.prompt.pos.X() * ibd.prompt.pos.X() + ibd.prompt.pos.Y() * ibd.prompt.pos.Y()) / 1.0e6, ibd.prompt.pos.Z() / 1000.0);
        h_rho_z_d->Fill((ibd.delayed.pos.X() * ibd.delayed.pos.X() + ibd.delayed.pos.Y() * ibd.delayed.pos.Y()) / 1.0e6, ibd.delayed.pos.Z() / 1000.0);
    }

    TCanvas* c_e_p = new TCanvas("c_e_p", "Prompt energy", 1000, 1000);
    c_e_p->cd();
    h_e_p->Draw();
    c_e_p->Update();

    TCanvas* c_e_d = new TCanvas("c_e_d", "Delayed energy", 1000, 1000);
    c_e_d->cd();
    h_e_d->Draw();
    c_e_d->Update();

    TCanvas* c_dt = new TCanvas("c_dt", "Prompt-Delayed time difference", 1000, 1000);
    c_dt->cd();
    h_dt->Draw();
    c_dt->Update();

    TCanvas* c_dr = new TCanvas("c_dr", "Prompt-Delayed distance", 1000, 1000);
    c_dr->cd();
    h_dr->Draw();
    c_dr->Update();

    TCanvas* c_rho_z_p = new TCanvas("c_rho_z_p", "Prompt vertex distribution", 1000, 1000);
    c_rho_z_p->cd();
    h_rho_z_p->Draw();
    c_rho_z_p->Update();

    TCanvas* c_rho_z_d = new TCanvas("c_rho_z_d", "Delayed vertex distribution", 1000, 1000);
    c_rho_z_d->cd();
    h_rho_z_d->Draw();
    c_rho_z_d->Update();
#endif
}