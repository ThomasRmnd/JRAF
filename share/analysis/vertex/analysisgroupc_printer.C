#include <algorithm>
#include <chrono>
#include <iostream>
#include <set>
#include <string>

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include "event.hpp"
#include "timestamp.hpp"
#include "vec3.hpp"

struct IBD {

    vertex prompt;
    vertex delayed;

};

struct Cosmo {

    vertex prompt;
    vertex delayed;
    double dlat_mu2p;
    double dlat_mu2d;
    double dt_mu2p;
    double dt_mu2d;

};

inline bool operator<(const IBD& lhs, const IBD& rhs) {
    return lhs.prompt.ts < rhs.prompt.ts;
}

inline bool operator<(const Cosmo& lhs, const Cosmo& rhs) {
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
        std::cout << "Prompt: " /* (" << posx_p << ", " << posy_p << ", " << posz_p << "), */ "E = " << e_p << ", Q = " << totq_p << ", Time = " << timestamp(sec_p, nsec_p) << '\n';
        std::cout << "Delayed: " /* (" << posx_d << ", " << posy_d << ", " << posz_d << "), */ "E = " << e_d << ", Q = " << totq_d << ", Time = " << timestamp(sec_d, nsec_d) << '\n';
    }

};

class MainAnalysis : public AnalysisBase {

public:

    MainAnalysis(const std::string& suffix) : AnalysisBase{"IBDAnalysis" + suffix} {}

    ~MainAnalysis() override = default;

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

    virtual bool selection() override = 0;

    void print() override {
        timestamp ts_p{sec_p, nsec_p};
        timestamp ts_d{sec_d, nsec_d};
        vec3 pos_p{posx_p, posy_p, posz_p};
        vec3 pos_d{posx_d, posy_d, posz_d};

        std::string method = "";
        double d_mu2p = std::numeric_limits<double>::infinity();
        double t_mu2p = std::numeric_limits<double>::infinity();
        double d_mu2d = std::numeric_limits<double>::infinity();
        double t_mu2d = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0ul; k < method_mu->size(); ++k) {
            timestamp ts_mu{sec_mu->operator[](k), nsec_mu->operator[](k)};
            vec3 pos_mu{posx_mu->operator[](k), posy_mu->operator[](k), posz_mu->operator[](k)};
            vec3 dir_mu{dirx_mu->operator[](k), diry_mu->operator[](k), dirz_mu->operator[](k)};
            double tmp_d_mu2p = mag(cross(dir_mu, pos_p - pos_mu));
            double tmp_t_mu2p = timestamp_to_double(ts_p - ts_mu);
            double tmp_d_mu2d = mag(cross(dir_mu, pos_d - pos_mu));
            double tmp_t_mu2d = timestamp_to_double(ts_d - ts_mu);
            if (tmp_d_mu2p < d_mu2p && 0.0 < tmp_t_mu2p && tmp_t_mu2p < 1.2) {
                method = method_mu->operator[](k);
                d_mu2p = tmp_d_mu2p;
                t_mu2p = tmp_t_mu2p;
                d_mu2d = tmp_d_mu2d;
                t_mu2d = tmp_t_mu2d;
            }
        }

        // std::cout << "Prompt: (" << posx_p << ", " << posy_p << ", " << posz_p << "), E = " << e_p << ", Q = " << totq_p << ", Time = " << TTimeStamp(sec_p, nsec_p) << '\n';
        // std::cout << "Delayed: (" << posx_d << ", " << posy_d << ", " << posz_d << "), E = " << e_d << ", Q = " << totq_d << ", Time = " << TTimeStamp(sec_d, nsec_d) << '\n';
        // std::cout << "Number of Neutron Veto associated: " << nb_neutron_veto << '\n';
        // std::cout << "Number of Multiplicity Veto associated: " << nb_multu_veto << '\n';
        // std::cout << "Muon (" << method << "): d_mu2p = " << d_mu2p << ", t_mu2p = " << t_mu2p << ", d_mu2d = " << d_mu2d << ", t_mu2d = " << t_mu2d << '\n';

        if (selection()) {
            std::cout << "Prompt: E = " << e_p << ", Q = " << totq_p << ", Time = " << timestamp{sec_p, nsec_p} << '\n';
            std::cout << "Delayed: E = " << e_d << ", Q = " << totq_d << ", Time = " << timestamp{sec_d, nsec_d} << '\n';
        }
    }

};

class IBDAnalysis : public MainAnalysis {

public:

    IBDAnalysis(const std::string& suffix) : MainAnalysis{suffix} {}

    ~IBDAnalysis() override = default;

    bool selection() override {
        if (stdt_p > 200.0 || stdt_d > 200.0) return false;
        timestamp ts_p{sec_p, nsec_p};
        timestamp ts_d{sec_d, nsec_d};
        vec3 pos_p{posx_p, posy_p, posz_p};
        vec3 pos_d{posx_d, posy_d, posz_d};
        if (mag(pos_p) > 16500.0) return false;
        if ((posz_p < -15500.0 || 15500 < posz_p) && std::sqrt(posx_p * posx_p + posy_p * posy_p) < 3000.0) return false;
        if (e_p < 0.7 || 12.0 < e_p) return false;
        if (e_d < 2.0 || 2.5 < e_d) return false;

        std::size_t nb_neutron_veto = 0ul;
        for (std::size_t k = 0ul; k < e_n->size(); ++k) {
            if (e_n->operator[](k) < 1.5 || 20.0 < e_n->operator[](k)) continue;
            timestamp ts_n{sec_n->operator[](k), nsec_n->operator[](k)};
            vec3 pos_n{posx_n->operator[](k), posy_n->operator[](k), posz_n->operator[](k)};
            // if (pos_n.Mag() > 17700.0) continue;
            if (mag(pos_p - pos_n) > 4000.0 || mag(pos_p - pos_n) > 4000.0) continue;
            if (ts_p < ts_n + timestamp{0, 20000} || ts_n + timestamp{0, 1200000000} < ts_p) continue;
            if (ts_d < ts_n + timestamp{0, 20000} || ts_n + timestamp{0, 1200000000} < ts_d) continue;
            ++nb_neutron_veto;
        }

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < e_mult->size(); ++k) {
            if (e_mult->operator[](k) < 2.0 || 12.0 < e_mult->operator[](k)) continue;
            timestamp ts_mult{sec_mult->operator[](k), nsec_mult->operator[](k)};
            vec3 pos_mult{posx_mult->operator[](k), posy_mult->operator[](k), posz_mult->operator[](k)};
            // if (pos_mult.Mag() > 17700.0) continue;
            // if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_mult < ts_p - timestamp{0, 1000000} || ts_d + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }

        return (nb_neutron_veto == 0ul && nb_multu_veto == 0ul);
    }

};

class CosmoAnalysis : public MainAnalysis {

public:

    CosmoAnalysis(const std::string& suffix, const timestamp& ts_low, const timestamp& ts_high, double radius) :
        MainAnalysis{suffix}, 
        m_ts_low{ts_low}, 
        m_ts_high{ts_high},
        m_radius{radius}
    {}

    ~CosmoAnalysis() override = default;

    bool selection() override {
        if (stdt_p > 200.0 || stdt_d > 200.0) return false;
        timestamp ts_p{sec_p, nsec_p};
        timestamp ts_d{sec_d, nsec_d};
        vec3 pos_p{posx_p, posy_p, posz_p};
        vec3 pos_d{posx_d, posy_d, posz_d};
        if (mag(pos_p) > 16500.0) return false;
        if ((posz_p < -15500.0 || 15500 < posz_p) && std::sqrt(posx_p * posx_p + posy_p * posy_p) < 3000.0) return false;
        if (e_p < 0.7 || 12.0 < e_p) return false;
        if (e_d < 2.0 || 2.5 < e_d) return false;

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < e_mult->size(); ++k) {
            if (e_mult->operator[](k) < 2.0 || 12.0 < e_mult->operator[](k)) continue;
            timestamp ts_mult{sec_mult->operator[](k), nsec_mult->operator[](k)};
            vec3 pos_mult{posx_mult->operator[](k), posy_mult->operator[](k), posz_mult->operator[](k)};
            // if (pos_mult.Mag() > 17700.0) continue;
            // if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_mult < ts_p - timestamp{0, 1000000} || ts_d + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }

        bool found = false;
        for (std::size_t k = 0ul; k < method_mu->size() && !found; ++k) {
            if (method_mu->operator[](k) != "CdWpTtChi2") continue;
            timestamp ts_mu{sec_mu->operator[](k), nsec_mu->operator[](k)};
            if (ts_p < ts_mu + m_ts_low || m_ts_high + ts_mu < ts_p) continue;
            if (ts_d < ts_mu + m_ts_low || m_ts_high + ts_mu < ts_d) continue;
            vec3 pos_mu{posx_mu->operator[](k), posy_mu->operator[](k), posz_mu->operator[](k)};
            vec3 dir_mu{dirx_mu->operator[](k), diry_mu->operator[](k), dirz_mu->operator[](k)};
            double d_mu2p = mag(cross(dir_mu, pos_p - pos_mu));
            double d_mu2d = mag(cross(dir_mu, pos_d - pos_mu));
            if (d_mu2p < m_radius && d_mu2d < m_radius) {
                found = true;
                m_dlat_mu2p = d_mu2p;
                m_dlat_mu2d = d_mu2d;
                m_dt_mu2p = timestamp_to_double(ts_p - ts_mu);
                m_dt_mu2d = timestamp_to_double(ts_d - ts_mu);
            }
        }

        return (nb_multu_veto == 0ul && found);
    }

    double dlat_p() const { return m_dlat_mu2p; }
    double dlat_d() const { return m_dlat_mu2d; }
    double dt_mu2p() const { return m_dt_mu2p; }
    double dt_mu2d() const { return m_dt_mu2d; }

private:

    timestamp m_ts_low;
    timestamp m_ts_high;
    double m_radius;

    double m_dlat_mu2p;
    double m_dlat_mu2d;
    double m_dt_mu2p;
    double m_dt_mu2d;

};


class CosmoRateWithNeutronAnalysis : public MainAnalysis {

public:

    CosmoRateWithNeutronAnalysis(const std::string& suffix) : MainAnalysis{suffix} {}

    ~CosmoRateWithNeutronAnalysis() override = default;

    bool selection() override {
        if (stdt_p > 200.0 || stdt_d > 200.0) return false;
        timestamp ts_p{sec_p, nsec_p};
        timestamp ts_d{sec_d, nsec_d};
        vec3 pos_p{posx_p, posy_p, posz_p};
        vec3 pos_d{posx_d, posy_d, posz_d};
        if (mag(pos_p) > 16500.0) return false;
        if ((posz_p < -15500.0 || 15500 < posz_p) && std::sqrt(posx_p * posx_p + posy_p * posy_p) < 3000.0) return false;
        if (e_p < 0.7 || 12.0 < e_p) return false;
        if (e_d < 2.0 || 2.5 < e_d) return false;

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < e_mult->size(); ++k) {
            if (e_mult->operator[](k) < 2.0 || 12.0 < e_mult->operator[](k)) continue;
            timestamp ts_mult{sec_mult->operator[](k), nsec_mult->operator[](k)};
            vec3 pos_mult{posx_mult->operator[](k), posy_mult->operator[](k), posz_mult->operator[](k)};
            // if (pos_mult.Mag() > 17700.0) continue;
            // if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_mult < ts_p - timestamp{0, 1000000} || ts_d + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }

        return (nb_multu_veto == 0ul);
    }

};

std::vector<IBD> get_all_ibd(const std::string& filename, AnalysisBase* analysis) {
    TChain* chain = analysis->retrieve(filename);
    if (!chain) return {};
    std::set<IBD> ibds_ordered;
    std::cout << "=== Analysis: " << analysis->name << " (Total Entries: " << chain->GetEntries() << ") ===\n";
    for (long k = 0; k < chain->GetEntries(); ++k) {
        chain->GetEntry(k);
        if (k % 1000 == 0) {
            std::cout << "\rProcessing Entry " << k << " / " << chain->GetEntries();
        }
        // analysis->print();
        if (!analysis->selection()) continue;
        IBD ibd;
        ibd.prompt.pos = vec3{analysis->posx_p, analysis->posy_p, analysis->posz_p};
        ibd.prompt.ts = timestamp{analysis->sec_p, analysis->nsec_p};
        ibd.prompt.e = analysis->e_p;
        ibd.prompt.q = analysis->totq_p;
        ibd.delayed.pos = vec3{analysis->posx_d, analysis->posy_d, analysis->posz_d};
        ibd.delayed.ts = timestamp{analysis->sec_d, analysis->nsec_d};
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

std::vector<Cosmo> get_all_cosmo(const std::string& filename, CosmoAnalysis* analysis) {
    TChain* chain = analysis->retrieve(filename);
    if (!chain) return {};
    std::set<Cosmo> cosmos_ordered;
    std::cout << "=== Analysis: Cosmo (Total Entries: " << chain->GetEntries() << ") ===\n";
    for (long k = 0; k < chain->GetEntries(); ++k) {
        chain->GetEntry(k);
        if (k % 1000 == 0) {
            std::cout << "\rProcessing Entry " << k << " / " << chain->GetEntries();
        }
        // analysis->print();
        if (!analysis->selection()) continue;
        Cosmo cosmo;
        cosmo.prompt.pos = vec3{analysis->posx_p, analysis->posy_p, analysis->posz_p};
        cosmo.prompt.ts = timestamp{analysis->sec_p, analysis->nsec_p};
        cosmo.prompt.e = analysis->e_p;
        cosmo.prompt.q = analysis->totq_p;
        cosmo.delayed.pos = vec3{analysis->posx_d, analysis->posy_d, analysis->posz_d};
        cosmo.delayed.ts = timestamp{analysis->sec_d, analysis->nsec_d};
        cosmo.delayed.e = analysis->e_d;
        cosmo.delayed.q = analysis->totq_d;
        cosmo.dlat_mu2p = analysis->dlat_p();
        cosmo.dlat_mu2d = analysis->dlat_d();
        cosmo.dt_mu2p = analysis->dt_mu2p();
        cosmo.dt_mu2d = analysis->dt_mu2d();
        cosmos_ordered.insert(cosmo);
    }

    std::vector<Cosmo> cosmos;
    cosmos.reserve(cosmos_ordered.size());
    for (std::set<Cosmo>::const_iterator it = cosmos_ordered.begin(); it != cosmos_ordered.end(); ++it) {
        cosmos.push_back(*it);
    }
    return cosmos;
}

void fit_and_plot_cosmo_rate_with_neutron(TH1D* h) {
    TCanvas* c = new TCanvas(Form("%s_canvas", h->GetName()), Form("%s Canvas", h->GetName()), 800, 600);
    c->cd();

    double constant_term = 0.0;
    for (int bin = h->GetXaxis()->FindBin(0.8); bin <= h->GetXaxis()->GetNbins(); ++bin) {
        constant_term += h->GetBinContent(bin);
    }
    constant_term /= (h->GetXaxis()->GetNbins() - h->GetXaxis()->FindBin(0.8) + 1);
    double exponential_term = h->GetMaximum() - constant_term;

    TF1* f = new TF1(Form("f_%s", h->GetName()), "[0] + [1] * exp(-x / [2])", 0.05, 1.2);
    f->SetParameter(0, constant_term);
    f->SetParameter(1, exponential_term);
    f->SetParameter(2, 180.0e-3);

    f->SetLineColor(kRed);
    f->SetLineWidth(3);

    h->Fit(f, "R");
    h->SetLineWidth(3);
    h->GetXaxis()->SetTitle("#Delta t_{#mu2p} (s)");
    h->GetXaxis()->CenterTitle(kTRUE);
    h->GetYaxis()->SetTitle("Entries");
    h->GetYaxis()->CenterTitle(kTRUE);
    h->GetYaxis()->SetTitleOffset(1.5);
    h->Draw();
    f->Draw("SAME");
    
    c->SetTickx();
    c->SetTicky();

    c->Update();

    double time_window = 1.2;
    double binning = time_window / 120.0;
    std::cout << "Fit Results for " << h->GetName() << ":\n";
    std::cout << "nIBD = " << f->GetParameter(0) * time_window / binning << " +/- " << f->GetParError(0) * time_window / binning << '\n';
    std::cout << "nLiHe = " << f->GetParameter(1) * f->GetParameter(2) * (1 - std::exp(-time_window / f->GetParameter(2))) / binning << '\n';    
}

void analyze_cosmo_rate_with_neutron(const std::string& filename, CosmoRateWithNeutronAnalysis* analysis) {
    TChain* chain = analysis->retrieve(filename);
    if (!chain) return;
    std::map<IBD, std::vector<double>> ibds_dt_mu2p;
    std::map<IBD, std::vector<double>> ibds_d_mu2p_cdwp;
    std::map<IBD, std::vector<double>> ibds_d_mu2p_tt;
    std::map<IBD, int> ibds_neutron_count;
    std::cout << "=== Analysis: CosmoRateWithNeutron (Total Entries: " << chain->GetEntries() << ") ===\n";
    for (long k = 0; k < chain->GetEntries(); ++k) {
        chain->GetEntry(k);
        if (k % 1000 == 0) {
            std::cout << "\rProcessing Entry " << k << " / " << chain->GetEntries();
        }
        // analysis->print();
        if (!analysis->selection()) continue;
        IBD ibd;
        ibd.prompt.pos = vec3{analysis->posx_p, analysis->posy_p, analysis->posz_p};
        ibd.prompt.ts = timestamp{analysis->sec_p, analysis->nsec_p};
        ibd.prompt.e = analysis->e_p;
        ibd.prompt.q = analysis->totq_p;
        ibd.delayed.pos = vec3{analysis->posx_d, analysis->posy_d, analysis->posz_d};
        ibd.delayed.ts = timestamp{analysis->sec_d, analysis->nsec_d};
        ibd.delayed.e = analysis->e_d;
        ibd.delayed.q = analysis->totq_d;

        std::vector<double> dt_mu2p_times;
        std::vector<double> d_mu2p_cdwp_values;
        std::vector<double> d_mu2p_tt_values;
        std::vector<std::pair<timestamp, std::vector<std::string>>> used_muon_times;
        int neutron_count = 0;
        for (std::size_t k = 0ul; k < analysis->method_mu->size(); ++k) {
            timestamp ts_mu{analysis->sec_mu->operator[](k), analysis->nsec_mu->operator[](k)};
            if (ibd.prompt.ts < ts_mu + timestamp{0, 5000000} || ts_mu + timestamp{0, 1200000000} < ibd.prompt.ts) continue;
            std::vector<std::pair<timestamp, std::vector<std::string>>>::iterator it = std::find_if(
                used_muon_times.begin(), used_muon_times.end(),
                [ts_mu](const std::pair<timestamp, std::vector<std::string>>& used_ts_mu) {
                    timestamp diff = ts_mu - used_ts_mu.first;
                    return (timestamp{0, -1000} < diff && diff < timestamp{0, 1000});
                }
            );
            if (it != used_muon_times.end()) {
                std::vector<std::string>::iterator method_it = std::find(
                    it->second.begin(), it->second.end(),
                    analysis->method_mu->operator[](k)
                );
                if (method_it != it->second.end()) continue;
                else {
                    it->second.push_back(analysis->method_mu->operator[](k));
                }
            }
            else {
                used_muon_times.push_back(std::make_pair(ts_mu, std::vector<std::string>{analysis->method_mu->operator[](k)}));
            }
            dt_mu2p_times.push_back(timestamp_to_double(ibd.prompt.ts - ts_mu));
            for (std::size_t i = 0ul; i < analysis->sec_n->size(); ++i) {
                if (analysis->e_n->operator[](i) < 2.0 || 2.5 < analysis->e_n->operator[](i)) continue;
                timestamp ts_n{analysis->sec_n->operator[](i), analysis->nsec_n->operator[](i)};
                if (ts_n < ts_mu + timestamp{0, 20000} || ts_mu + timestamp{0, 2000000} < ts_n) continue;
                ++neutron_count;
            }
            vec3 pos_mu{analysis->posx_mu->operator[](k), analysis->posy_mu->operator[](k), analysis->posz_mu->operator[](k)};
            vec3 dir_mu{analysis->dirx_mu->operator[](k), analysis->diry_mu->operator[](k), analysis->dirz_mu->operator[](k)};
            double d_mu2p = mag(cross(dir_mu, ibd.prompt.pos - pos_mu));
            if (analysis->method_mu->operator[](k) == "CdWpTtChi2") {
                d_mu2p_cdwp_values.push_back(d_mu2p);
            }
            if (analysis->method_mu->operator[](k) == "Tt") {
                d_mu2p_tt_values.push_back(d_mu2p);
            }
        }
        ibds_dt_mu2p[ibd] = dt_mu2p_times;
        ibds_d_mu2p_cdwp[ibd] = d_mu2p_cdwp_values;
        ibds_d_mu2p_tt[ibd] = d_mu2p_tt_values;
        ibds_neutron_count[ibd] = neutron_count;
    }

    TH1D* h_cosmo_rate_with_neutron = new TH1D("h_cosmo_rate_with_neutron", "Cosmo Rate With Neutron", 120, 0.0, 1.2);
    TH1D* h_cosmo_rate_with_no_neutron = new TH1D("h_cosmo_rate_with_no_neutron", "Cosmo Rate With No Neutron", 120, 0.0, 1.2);
    TH1D* h_cosmo_rate_with_at_least_1_neutron = new TH1D("h_cosmo_rate_with_at_least_1_neutron", "Cosmo Rate With At Least 1 Neutron", 120, 0.0, 1.2);
    TH1D* h_cosmo_rate_with_at_least_2_neutron = new TH1D("h_cosmo_rate_with_at_least_2_neutron", "Cosmo Rate With At Least 2 Neutron", 120, 0.0, 1.2);
    TH1D* h_cosmo_rate_with_at_least_3_neutron = new TH1D("h_cosmo_rate_with_at_least_3_neutron", "Cosmo Rate With At Least 3 Neutron", 120, 0.0, 1.2);
    for (const std::pair<IBD, std::vector<double>>& ibd_dt_mu2p : ibds_dt_mu2p) {
        for (double dt_mu2p : ibd_dt_mu2p.second) {
            h_cosmo_rate_with_neutron->Fill(dt_mu2p);
        }
        int neutron_count = ibds_neutron_count[ibd_dt_mu2p.first];
        if (neutron_count == 0) {
            for (double dt_mu2p : ibd_dt_mu2p.second) {
                h_cosmo_rate_with_no_neutron->Fill(dt_mu2p);
            }
        }
        if (neutron_count >= 1) {
            for (double dt_mu2p : ibd_dt_mu2p.second) {
                h_cosmo_rate_with_at_least_1_neutron->Fill(dt_mu2p);
            }
        }
        if (neutron_count >= 2) {
            for (double dt_mu2p : ibd_dt_mu2p.second) {
                h_cosmo_rate_with_at_least_2_neutron->Fill(dt_mu2p);
            }
        }
        if (neutron_count >= 3) {
            for (double dt_mu2p : ibd_dt_mu2p.second) {
                h_cosmo_rate_with_at_least_3_neutron->Fill(dt_mu2p);
            }
        }
    }

    fit_and_plot_cosmo_rate_with_neutron(h_cosmo_rate_with_neutron);
    fit_and_plot_cosmo_rate_with_neutron(h_cosmo_rate_with_no_neutron);
    fit_and_plot_cosmo_rate_with_neutron(h_cosmo_rate_with_at_least_1_neutron);
    fit_and_plot_cosmo_rate_with_neutron(h_cosmo_rate_with_at_least_2_neutron);
    fit_and_plot_cosmo_rate_with_neutron(h_cosmo_rate_with_at_least_3_neutron);

    TH2D* h_d_mu2p_cdwp_vs_dt_mu2p = new TH2D("h_d_mu2p_cdwp_vs_dt_mu2p", "Cosmo time vs distance", 120, 0.0, 1.5, 100, 0.0, 40000.0);
    TH2D* h_d_mu2p_tt_vs_dt_mu2p = new TH2D("h_d_mu2p_tt_vs_dt_mu2p", "Cosmo time vs distance", 120, 0.0, 1.5, 100, 0.0, 40000.0);
    for (const std::pair<IBD, std::vector<double>>& ibd_dt_mu2p : ibds_dt_mu2p) {
        const IBD& ibd = ibd_dt_mu2p.first;
        const std::vector<double>& dt_mu2p_times = ibd_dt_mu2p.second;
        const std::vector<double>& d_mu2p_cdwp_values = ibds_d_mu2p_cdwp[ibd];
        const std::vector<double>& d_mu2p_tt_values = ibds_d_mu2p_tt[ibd];
        for (std::size_t k = 0ul; k < dt_mu2p_times.size(); ++k) {
            h_d_mu2p_cdwp_vs_dt_mu2p->Fill(dt_mu2p_times[k], d_mu2p_cdwp_values[k]);
            h_d_mu2p_tt_vs_dt_mu2p->Fill(dt_mu2p_times[k], d_mu2p_tt_values[k]);
        }
    }

    TCanvas* c_d_mu2p_cdwp_vs_dt_mu2p = new TCanvas("c_d_mu2p_cdwp_vs_dt_mu2p", "d_mu2p_cdwp vs dt_mu2p", 1000, 1000);
    c_d_mu2p_cdwp_vs_dt_mu2p->cd();
    h_d_mu2p_cdwp_vs_dt_mu2p->GetXaxis()->SetTitle("#Delta t_{#mu2p} (s)");
    h_d_mu2p_cdwp_vs_dt_mu2p->GetXaxis()->CenterTitle(kTRUE);
    h_d_mu2p_cdwp_vs_dt_mu2p->GetYaxis()->SetTitle("d_{#mu2p} (cm)");
    h_d_mu2p_cdwp_vs_dt_mu2p->GetYaxis()->CenterTitle(kTRUE);
    h_d_mu2p_cdwp_vs_dt_mu2p->GetYaxis()->SetTitleOffset(1.5);
    h_d_mu2p_cdwp_vs_dt_mu2p->Draw("COLZ");
    c_d_mu2p_cdwp_vs_dt_mu2p->SetTickx();
    c_d_mu2p_cdwp_vs_dt_mu2p->SetTicky();
    c_d_mu2p_cdwp_vs_dt_mu2p->Update();

    TCanvas* c_d_mu2p_tt_vs_dt_mu2p = new TCanvas("c_d_mu2p_tt_vs_dt_mu2p", "d_mu2p_tt vs dt_mu2p", 1000, 1000);
    c_d_mu2p_tt_vs_dt_mu2p->cd();
    h_d_mu2p_tt_vs_dt_mu2p->GetXaxis()->SetTitle("#Delta t_{#mu2p} (s)");
    h_d_mu2p_tt_vs_dt_mu2p->GetXaxis()->CenterTitle(kTRUE);
    h_d_mu2p_tt_vs_dt_mu2p->GetYaxis()->SetTitle("d_{#mu2p} (cm)");
    h_d_mu2p_tt_vs_dt_mu2p->GetYaxis()->CenterTitle(kTRUE);
    h_d_mu2p_tt_vs_dt_mu2p->GetYaxis()->SetTitleOffset(1.5);
    h_d_mu2p_tt_vs_dt_mu2p->Draw("COLZ");
    c_d_mu2p_tt_vs_dt_mu2p->SetTickx();
    c_d_mu2p_tt_vs_dt_mu2p->SetTicky();
    c_d_mu2p_tt_vs_dt_mu2p->Update();
}

void print_all_entries(const std::string& filename, AnalysisBase* analysis) {
    TChain* chain = analysis->retrieve(filename);
    if (!chain) return;
    std::set<IBD> ibds;
    std::cout << "=== Analysis: " << analysis->name << " (Total Entries: " << chain->GetEntries() << ") ===\n";
    for (long k = 0; k < chain->GetEntries(); ++k) {
        chain->GetEntry(k);
        if (k % 1000 == 0) {
            std::cout << "\rProcessing Entry " << k << " / " << chain->GetEntries();
        }
        // analysis->print();
        if (!analysis->selection()) continue;
        IBD ibd;
        ibd.prompt.pos = vec3{analysis->posx_p, analysis->posy_p, analysis->posz_p};
        ibd.prompt.ts = timestamp{analysis->sec_p, analysis->nsec_p};
        ibd.prompt.e = analysis->e_p;
        ibd.prompt.q = analysis->totq_p;
        ibd.delayed.pos = vec3{analysis->posx_d, analysis->posy_d, analysis->posz_d};
        ibd.delayed.ts = timestamp{analysis->sec_d, analysis->nsec_d};
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

    timestamp ts_p;
    double e_p;
    double totq_p;

    timestamp ts_d;
    double e_d;
    double totq_d;

};

bool operator<(const VanessaIBD& lhs, const VanessaIBD& rhs) {
    return lhs.ts_p < rhs.ts_p;
}

std::vector<VanessaIBD> analyze_vanessa_result() {
    TFile* file = TFile::Open("/sps/juno/jdeandre/rtraw_ThomasRaymond/analysis/ibd/summary/ReProd25C/IBD_all_reprodC.root", "READ");
    TTree* tree = file->Get<TTree>("events");
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
        ibd.ts_p = timestamp{sec_p, nsec_p};
        ibd.e_p = e_p;
        ibd.totq_p = totq_p;
        ibd.ts_d = timestamp{sec_d, nsec_d};
        ibd.e_d = e_d;
        ibd.totq_d = totq_d;
        ibds.insert(ibd);
    }
    std::vector<VanessaIBD> ibds_vector;
    ibds_vector.reserve(ibds.size());
    for (std::set<VanessaIBD>::const_iterator it = ibds.begin(); it != ibds.end(); ++it) {
        ibds_vector.push_back(*it);
    }
    return ibds_vector;
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
        ibd.ts_p = timestamp{sec_p, nsec_p};
        ibd.e_p = e_p;
        ibd.totq_p = totq_p;
        ibd.ts_d = timestamp{sec_d, nsec_d};
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
        if (k % 1000 == 0) {
            std::cout << "\rProcessing Entry " << k << " / " << chain->GetEntries();
        }
        std::set<VanessaIBD>::const_iterator it = std::find_if(
            vanessa_ibds.begin(),
            vanessa_ibds.end(),
            [&](const VanessaIBD& vanessa_ibd) {
                return (
                    vanessa_ibd.ts_p.sec == analysis->sec_p &&
                    vanessa_ibd.ts_d.sec == analysis->sec_d &&
                    vanessa_ibd.e_p == analysis->e_p &&
                    vanessa_ibd.e_d == analysis->e_d
                );
            }
        );
        bool is_only_in_vanessa = (it != vanessa_ibds.end() && !analysis->selection());
        bool is_only_in_analysis = (it == vanessa_ibds.end() && analysis->selection());
        if (!is_only_in_vanessa && !is_only_in_analysis) continue;

        timestamp ts_p{analysis->sec_p, analysis->nsec_p};
        timestamp ts_d{analysis->sec_d, analysis->nsec_d};
        vec3 pos_p{analysis->posx_p, analysis->posy_p, analysis->posz_p};
        vec3 pos_d{analysis->posx_d, analysis->posy_d, analysis->posz_d};

        std::size_t nb_neutron_veto = 0ul;
        for (std::size_t k = 0ul; k < analysis->e_n->size(); ++k) {
            if (analysis->e_n->operator[](k) < 1.5 || 20.0 < analysis->e_n->operator[](k)) continue;
            timestamp ts_n{analysis->sec_n->operator[](k), analysis->nsec_n->operator[](k)};
            vec3 pos_n{analysis->posx_n->operator[](k), analysis->posy_n->operator[](k), analysis->posz_n->operator[](k)};
            // if (pos_n.Mag() > 17700.0) continue;
            if (mag(pos_p - pos_n) > 4000.0 || mag(pos_p - pos_n) > 4000.0) continue;
            if (ts_p - ts_n < timestamp{0, 20000} || timestamp{0, 1200000000} < ts_p - ts_n) continue;
            if (ts_d - ts_n < timestamp{0, 20000} || timestamp{0, 1200000000} < ts_d - ts_n) continue;
            ++nb_neutron_veto;
        }

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < analysis->e_mult->size(); ++k) {
            if (analysis->e_mult->operator[](k) < 2.0 || 12.0 < analysis->e_mult->operator[](k)) continue;
            timestamp ts_mult{analysis->sec_mult->operator[](k), analysis->nsec_mult->operator[](k)};
            vec3 pos_mult{analysis->posx_mult->operator[](k), analysis->posy_mult->operator[](k), analysis->posz_mult->operator[](k)};
            // if (pos_mult.Mag() > 17700.0) continue;
            // if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_mult < ts_p - timestamp{0, 1000000} || ts_d + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }

        std::string method = "";
        double d_mu2p = std::numeric_limits<double>::infinity();
        double t_mu2p = std::numeric_limits<double>::infinity();
        double d_mu2d = std::numeric_limits<double>::infinity();
        double t_mu2d = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0ul; k < analysis->method_mu->size(); ++k) {
            timestamp ts_mu{analysis->sec_mu->operator[](k), analysis->nsec_mu->operator[](k)};
            vec3 pos_mu{analysis->posx_mu->operator[](k), analysis->posy_mu->operator[](k), analysis->posz_mu->operator[](k)};
            vec3 dir_mu{analysis->dirx_mu->operator[](k), analysis->diry_mu->operator[](k), analysis->dirz_mu->operator[](k)};
            double tmp_d_mu2p = mag(cross(dir_mu, pos_p - pos_mu));
            double tmp_t_mu2p = timestamp_to_double(ts_p - ts_mu);
            double tmp_d_mu2d = mag(cross(dir_mu, pos_d - pos_mu));
            double tmp_t_mu2d = timestamp_to_double(ts_d - ts_mu);
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
        std::cout << "Prompt: (" << analysis->posx_p << ", " << analysis->posy_p << ", " << analysis->posz_p << "), E = " << analysis->e_p << ", Q = " << analysis->totq_p << ", Time = " << timestamp{analysis->sec_p, analysis->nsec_p} << '\n';
        std::cout << "Delayed: (" << analysis->posx_d << ", " << analysis->posy_d << ", " << analysis->posz_d << "), E = " << analysis->e_d << ", Q = " << analysis->totq_d << ", Time = " << timestamp{analysis->sec_d, analysis->nsec_d} << '\n';
        std::cout << "Number of Neutron Veto associated: " << nb_neutron_veto << '\n';
        std::cout << "Number of Multiplicity Veto associated: " << nb_multu_veto << '\n';

        // print neutrons
        for (std::size_t k = 0ul; k < analysis->e_n->size(); ++k) {
            timestamp ts_n{analysis->sec_n->operator[](k), analysis->nsec_n->operator[](k)};
            vec3 pos_n{analysis->posx_n->operator[](k), analysis->posy_n->operator[](k), analysis->posz_n->operator[](k)};
            double e_n = analysis->e_n->operator[](k);
            double totq_n = analysis->totq_n->operator[](k);
            std::cout << "  Neutron: " << pos_n << ", E = " << e_n << ", Q = " << totq_n << ", Time = " << ts_n << '\n';
        }

        // print mults
        for (std::size_t k = 0ul; k < analysis->e_mult->size(); ++k) {
            timestamp ts_mult{analysis->sec_mult->operator[](k), analysis->nsec_mult->operator[](k)};
            vec3 pos_mult{analysis->posx_mult->operator[](k), analysis->posy_mult->operator[](k), analysis->posz_mult->operator[](k)};
            double e_mult = analysis->e_mult->operator[](k);
            double totq_mult = analysis->totq_mult->operator[](k);
            std::cout << "  Mult: " << pos_mult << ", E = " << e_mult << ", Q = " << totq_mult << ", Time = " << ts_mult << '\n';
        }

        // print muons
        for (std::size_t k = 0ul; k < analysis->method_mu->size(); ++k) {
            timestamp ts_mu{analysis->sec_mu->operator[](k), analysis->nsec_mu->operator[](k)};
            vec3 pos_mu{analysis->posx_mu->operator[](k), analysis->posy_mu->operator[](k), analysis->posz_mu->operator[](k)};
            vec3 dir_mu{analysis->dirx_mu->operator[](k), analysis->diry_mu->operator[](k), analysis->dirz_mu->operator[](k)};
            double totq_mu = analysis->totq_mu->operator[](k);
            double quality_mu = analysis->quality_mu->operator[](k);
            std::cout << "  Muon (" << analysis->method_mu->operator[](k) << "): " << pos_mu << ", Dir = " << dir_mu << ", Q = " << totq_mu << ", Quality = " << quality_mu << ", Time = " << ts_mu << '\n';
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

std::vector<double> generate_segment_boundaries(double start, double stop, int num_bins) {
    if (num_bins <= 0) return {};
    int num_points = num_bins + 1;
    double expected_width = (stop - start) / num_bins;
    
    std::vector<double> segment;
    segment.reserve(num_points);
    segment.push_back(start); 

    for (int i = 1; i < num_points; ++i) {
        double boundary = start + i * expected_width;
        
        if (i == num_points - 1) {
             segment.push_back(stop);
        } else {
             segment.push_back(boundary);
        }
    }
    return segment;
}

std::vector<double> create_custom_e_p_bins() {

    // double s1_start = 0.8;
    // double s2_start = 0.94;
    // double s3_start = 7.44;
    // double s4_start = 7.8;
    // double s5_start = 8.2;
    // double stop = 12.0;

    // int s1_bins = 1;
    // int s2_bins = 325;
    // int s3_bins = 9;
    // int s4_bins = 4;
    // int s5_bins = 1;
    // int tot_bins = s1_bins + s2_bins + s3_bins + s4_bins + s5_bins;

    double edges[] = {0.7, 1.0, 6.6, 7.4, 7.7, 8.1, 8.6, 9.4, 12.0};
    int    bins[]  = {  1,  56,   4,   1,   1,   1,   1,   1};
    
    std::vector<double> e_p_bins;

    for (std::size_t k = 0ul; k < 8ul; ++k) {
        double start = edges[k];
        double stop = edges[k + 1];
        int nbins = bins[k];
        std::vector<double> segment = generate_segment_boundaries(start, stop, nbins);
        if (k == 0ul) {
            e_p_bins.insert(e_p_bins.end(), segment.begin(), segment.end());
        } else {
            e_p_bins.insert(e_p_bins.end(), segment.begin() + 1, segment.end());
        }
    }
    
    return e_p_bins;
}

#undef PRINT_ALL_ENTRIES
#undef COMPARE_WITH_VANESSA
#define GET_ALL_IBD

void analysisgroupc_printer(const std::string& filename, const std::string& suffix) {

    CosmoRateWithNeutronAnalysis cosmo_rate_with_neutron_analysis(suffix);
    analyze_cosmo_rate_with_neutron(filename, &cosmo_rate_with_neutron_analysis);

    CosmoAnalysis cosmo_before_analysis(suffix, timestamp{0, -1200000000}, timestamp{0, -5000000}, 3000.0);
    CosmoAnalysis cosmo_after_analysis(suffix, timestamp{0, 5000000}, timestamp{0, 1200000000}, 3000.0);
    std::vector<Cosmo> cosmos_before = get_all_cosmo(filename, &cosmo_before_analysis);
    std::vector<Cosmo> cosmos_after = get_all_cosmo(filename, &cosmo_after_analysis);

    IBDAnalysis ibd_analysis(suffix);
#ifdef PRINT_ALL_ENTRIES
    print_all_entries(filename, &ibd_analysis);
#endif 
#ifdef COMPARE_WITH_VANESSA
    compare_with_vanessa(filename, &ibd_analysis);
#endif
#ifdef GET_ALL_IBD
    std::vector<IBD> ibds = get_all_ibd(filename, &ibd_analysis);
    std::vector<VanessaIBD> vanessa_ibds = analyze_vanessa_result();


    // double e_p_min = 0.7;
    // double e_p_max = 12.0;
    // double e_p_width = 0.20;
    // int e_p_nbin = std::round((e_p_max - e_p_min) / e_p_width) + 1;
    // std::vector<double> e_p_bins = linspace_cpp(e_p_min, e_p_max, e_p_nbin);
    std::vector<double> e_p_bins = create_custom_e_p_bins();
    TH1D* h_e_p = new TH1D("h_e_p", "Prompt energy", e_p_bins.size() - 1, e_p_bins.data());
    TH1D* h_e_p_vanessa = new TH1D("h_e_p_vanessa", "Prompt energy (Vanessa)", e_p_bins.size() - 1, e_p_bins.data());
    TH1D* h_e_p_cosmo_before = new TH1D("h_e_p_cosmo_before", "Prompt energy (Cosmo Before)", e_p_bins.size() - 1, e_p_bins.data());
    TH1D* h_e_p_cosmo_after = new TH1D("h_e_p_cosmo_after", "Prompt energy (Cosmo After)", e_p_bins.size() - 1, e_p_bins.data());
    TH1D* h_e_p_cosmo_diff = new TH1D("h_e_p_cosmo_diff", "Prompt energy (Cosmo After - Cosmo Before)", e_p_bins.size() - 1, e_p_bins.data());

    double e_d_min = 2.0;
    double e_d_max = 2.5;
    double e_d_width = 0.02;
    int e_d_nbin = std::round((e_d_max - e_d_min) / e_d_width) + 1;
    std::vector<double> e_d_bins = linspace_cpp(e_d_min, e_d_max, e_d_nbin);
    TH1D* h_e_d = new TH1D("h_e_d", "Delayed energy", e_d_bins.size() - 1, e_d_bins.data());
    TH1D* h_e_d_vanessa = new TH1D("h_e_d_vanessa", "Delayed energy (Vanessa)", e_d_bins.size() - 1, e_d_bins.data());
    TH1D* h_e_d_cosmo_before = new TH1D("h_e_d_cosmo_before", "Delayed energy (Cosmo Before)", e_d_bins.size() - 1, e_d_bins.data());
    TH1D* h_e_d_cosmo_after = new TH1D("h_e_d_cosmo_after", "Delayed energy (Cosmo After)", e_d_bins.size() - 1, e_d_bins.data());

    double e_dt_min = 0.0;
    double e_dt_max = 1.0;
    double e_dt_width = 0.025;
    int e_dt_nbin = std::round((e_dt_max - e_dt_min) / e_dt_width) + 1;
    std::vector<double> e_dt_bins = linspace_cpp(e_dt_min, e_dt_max, e_dt_nbin);
    TH1D* h_dt = new TH1D("h_dt", "Prompt-Delayed time difference", e_dt_bins.size() - 1, e_dt_bins.data());
    TH1D* h_dt_vanessa = new TH1D("h_dt_vanessa", "Prompt-Delayed time difference (Vanessa)", e_dt_bins.size() - 1, e_dt_bins.data());
    TH1D* h_dt_cosmo_before = new TH1D("h_dt_cosmo_before", "Prompt-Delayed time difference (Cosmo Before)", e_dt_bins.size() - 1, e_dt_bins.data());
    TH1D* h_dt_cosmo_after = new TH1D("h_dt_cosmo_after", "Prompt-Delayed time difference (Cosmo After)", e_dt_bins.size() - 1, e_dt_bins.data());

    double e_dr_min = 0.0;
    double e_dr_max = 1.5;
    double e_dr_width = 0.05;
    int e_dr_nbin = std::round((e_dr_max - e_dr_min) / e_dr_width) + 1;
    std::vector<double> e_dr_bins = linspace_cpp(e_dr_min, e_dr_max, e_dr_nbin);
    TH1D* h_dr = new TH1D("h_dr", "Prompt-Delayed distance", e_dr_bins.size() - 1, e_dr_bins.data());
    TH1D* h_dr_cosmo_before = new TH1D("h_dr_cosmo_before", "Prompt-Delayed distance (Cosmo Before)", e_dr_bins.size() - 1, e_dr_bins.data());
    TH1D* h_dr_cosmo_after = new TH1D("h_dr_cosmo_after", "Prompt-Delayed distance (Cosmo After)", e_dr_bins.size() - 1, e_dr_bins.data());

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
    TH2D* h_rho_z_p_cosmo_before = new TH2D("h_rho_z_p_cosmo_before", "Prompt vertex distribution (Cosmo Before)", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
    TH2D* h_rho_z_d_cosmo_before = new TH2D("h_rho_z_d_cosmo_before", "Delayed vertex distribution (Cosmo Before)", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
    TH2D* h_rho_z_p_cosmo_after = new TH2D("h_rho_z_p_cosmo_after", "Prompt vertex distribution (Cosmo After)", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());
    TH2D* h_rho_z_d_cosmo_after = new TH2D("h_rho_z_d_cosmo_after", "Delayed vertex distribution (Cosmo After)", rho_bins.size() - 1, rho_bins.data(), z_bins.size() - 1, z_bins.data());

    for (const IBD& ibd : ibds) {
        h_e_p->Fill(ibd.prompt.e);
        h_e_d->Fill(ibd.delayed.e);
        h_dt->Fill(timestamp_to_double(ibd.delayed.ts - ibd.prompt.ts) * 1000.0);
        h_dr->Fill(mag(ibd.delayed.pos - ibd.prompt.pos) / 1000.0);
        h_rho_z_p->Fill((ibd.prompt.pos.x * ibd.prompt.pos.x + ibd.prompt.pos.y * ibd.prompt.pos.y) / 1.0e6, ibd.prompt.pos.z / 1000.0);
        h_rho_z_d->Fill((ibd.delayed.pos.x * ibd.delayed.pos.x + ibd.delayed.pos.y * ibd.delayed.pos.y) / 1.0e6, ibd.delayed.pos.z / 1000.0);
    }

    for (const Cosmo& cosmo : cosmos_before) {
        h_e_p_cosmo_before->Fill(cosmo.prompt.e);
        h_e_d_cosmo_before->Fill(cosmo.delayed.e);
        h_dt_cosmo_before->Fill(timestamp_to_double(cosmo.delayed.ts - cosmo.prompt.ts) * 1000.0);
        h_dr_cosmo_before->Fill(mag(cosmo.delayed.pos - cosmo.prompt.pos) / 1000.0);
        h_rho_z_p_cosmo_before->Fill((cosmo.prompt.pos.x * cosmo.prompt.pos.x + cosmo.prompt.pos.y * cosmo.prompt.pos.y) / 1.0e6, cosmo.prompt.pos.z / 1000.0);
        h_rho_z_d_cosmo_before->Fill((cosmo.delayed.pos.x * cosmo.delayed.pos.x + cosmo.delayed.pos.y * cosmo.delayed.pos.y) / 1.0e6, cosmo.delayed.pos.z / 1000.0);
    }

    for (const Cosmo& cosmo : cosmos_after) {
        h_e_p_cosmo_after->Fill(cosmo.prompt.e);
        h_e_d_cosmo_after->Fill(cosmo.delayed.e);
        h_dt_cosmo_after->Fill(timestamp_to_double(cosmo.delayed.ts - cosmo.prompt.ts) * 1000.0);
        h_dr_cosmo_after->Fill(mag(cosmo.delayed.pos - cosmo.prompt.pos) / 1000.0);
        h_rho_z_p_cosmo_after->Fill((cosmo.prompt.pos.x * cosmo.prompt.pos.x + cosmo.prompt.pos.y * cosmo.prompt.pos.y) / 1.0e6, cosmo.prompt.pos.z / 1000.0);
        h_rho_z_d_cosmo_after->Fill((cosmo.delayed.pos.x * cosmo.delayed.pos.x + cosmo.delayed.pos.y * cosmo.delayed.pos.y) / 1.0e6, cosmo.delayed.pos.z / 1000.0);
    }

    h_e_p_cosmo_diff->Add(h_e_p_cosmo_after, h_e_p_cosmo_before, 1.0, -1.0);

    for (const VanessaIBD& ibd : vanessa_ibds) {
        h_e_p_vanessa->Fill(ibd.e_p);
        h_e_d_vanessa->Fill(ibd.e_d);
        h_dt_vanessa->Fill(timestamp_to_double(ibd.ts_d - ibd.ts_p) * 1000.0);
        // h_dr_vanessa->Fill((ibd.delayed.pos - ibd.prompt.pos).Mag() / 1000.0);
        // h_rho_z_p_vanessa->Fill((ibd.prompt.pos.X() * ibd.prompt.pos.X() + ibd.prompt.pos.Y() * ibd.prompt.pos.Y()) / 1.0e6, ibd.prompt.pos.Z() / 1000.0);
        // h_rho_z_d_vanessa->Fill((ibd.delayed.pos.X() * ibd.delayed.pos.X() + ibd.delayed.pos.Y() * ibd.delayed.pos.Y()) / 1.0e6, ibd.delayed.pos.Z() / 1000.0);
    }

    // ============================================================================================
    // Prompt energy
    // ============================================================================================

    TCanvas* c_e_p = new TCanvas("c_e_p", "Prompt energy", 1000, 1000);
    c_e_p->cd();

    h_e_p->SetLineWidth(3);
    h_e_p->SetLineStyle(kSolid);
    h_e_p->SetLineColorAlpha(kBlue, 1.0);

    h_e_p_vanessa->SetLineWidth(3);
    h_e_p_vanessa->SetLineStyle(kSolid);
    h_e_p_vanessa->SetLineColorAlpha(kRed, 1.0);

    h_e_p->Draw();
    h_e_p_vanessa->Draw("SAME");

    c_e_p->Update();

    // ============================================================================================
    // Delayed energy
    // ============================================================================================

    TCanvas* c_e_d = new TCanvas("c_e_d", "Delayed energy", 1000, 1000);
    c_e_d->cd();

    h_e_d->SetLineWidth(3);
    h_e_d->SetLineStyle(kSolid);
    h_e_d->SetLineColorAlpha(kBlue, 1.0);

    h_e_d_vanessa->SetLineWidth(3);
    h_e_d_vanessa->SetLineStyle(kSolid);
    h_e_d_vanessa->SetLineColorAlpha(kRed, 1.0);

    h_e_d->Draw();
    h_e_d_vanessa->Draw("SAME");

    c_e_d->Update();

    // ============================================================================================
    // Prompt-Delayed time difference
    // ============================================================================================

    TCanvas* c_dt = new TCanvas("c_dt", "Prompt-Delayed time difference", 1000, 1000);
    c_dt->cd();

    h_dt->SetLineWidth(3);
    h_dt->SetLineStyle(kSolid);
    h_dt->SetLineColorAlpha(kBlue, 1.0);

    h_dt_vanessa->SetLineWidth(3);
    h_dt_vanessa->SetLineStyle(kSolid);
    h_dt_vanessa->SetLineColorAlpha(kRed, 1.0);

    h_dt->Draw();
    h_dt_vanessa->Draw("SAME");
    
    c_dt->Update();

    // ============================================================================================
    // Prompt-Delayed distance
    // ============================================================================================

    TCanvas* c_dr = new TCanvas("c_dr", "Prompt-Delayed distance", 1000, 1000);
    c_dr->cd();
    h_dr->Draw();
    c_dr->Update();

    // ============================================================================================
    // Prompt vertex position
    // ============================================================================================

    TCanvas* c_rho_z_p = new TCanvas("c_rho_z_p", "Prompt vertex distribution", 1000, 1000);
    c_rho_z_p->cd();
    h_rho_z_p->Draw();
    c_rho_z_p->Update();

    // ============================================================================================
    // Delayed vertex position
    // ============================================================================================

    TCanvas* c_rho_z_d = new TCanvas("c_rho_z_d", "Delayed vertex distribution", 1000, 1000);
    c_rho_z_d->cd();
    h_rho_z_d->Draw();
    c_rho_z_d->Update();
#endif

    // ============================================================================================
    // Cosmo before - Prompt energy
    // ============================================================================================

    TCanvas* c_e_p_cosmo_before = new TCanvas("c_e_p_cosmo_before", "Prompt energy (Cosmo before)", 1000, 1000);
    c_e_p_cosmo_before->cd();

    h_e_p_cosmo_before->SetLineWidth(3);
    h_e_p_cosmo_before->SetLineStyle(kSolid);
    h_e_p_cosmo_before->SetLineColorAlpha(kBlue, 1.0);

    h_e_p_cosmo_before->Draw();

    c_e_p_cosmo_before->Update();

    // ============================================================================================
    // Cosmo before - Delayed energy
    // ============================================================================================

    TCanvas* c_e_d_cosmo_before = new TCanvas("c_e_d_cosmo_before", "Delayed energy (Cosmo before)", 1000, 1000);
    c_e_d_cosmo_before->cd();

    h_e_d_cosmo_before->SetLineWidth(3);
    h_e_d_cosmo_before->SetLineStyle(kSolid);
    h_e_d_cosmo_before->SetLineColorAlpha(kBlue, 1.0);

    h_e_d_cosmo_before->Draw();

    c_e_d_cosmo_before->Update();

    // ============================================================================================
    // Cosmo before - Prompt-Delayed time difference
    // ============================================================================================

    TCanvas* c_dt_cosmo_before = new TCanvas("c_dt_cosmo_before", "Prompt-Delayed time difference (Cosmo before)", 1000, 1000);
    c_dt_cosmo_before->cd();
    
    h_dt_cosmo_before->SetLineWidth(3);
    h_dt_cosmo_before->SetLineStyle(kSolid);
    h_dt_cosmo_before->SetLineColorAlpha(kBlue, 1.0);

    h_dt_cosmo_before->Draw();

    c_dt_cosmo_before->Update();

    // ============================================================================================
    // Cosmo before - Prompt-Delayed distance
    // ============================================================================================

    TCanvas* c_dr_cosmo_before = new TCanvas("c_dr_cosmo_before", "Prompt-Delayed distance (Cosmo before)", 1000, 1000);
    c_dr_cosmo_before->cd();

    h_dr_cosmo_before->Draw();

    c_dr_cosmo_before->Update();

    // ============================================================================================
    // Cosmo before - Prompt vertex position
    // ============================================================================================

    TCanvas* c_rho_z_p_cosmo_before = new TCanvas("c_rho_z_p_cosmo_before", "Prompt vertex distribution (Cosmo before)", 1000, 1000);
    c_rho_z_p_cosmo_before->cd();

    h_rho_z_p_cosmo_before->Draw();

    c_rho_z_p_cosmo_before->Update();

    // ============================================================================================
    // Cosmo before - Delayed vertex position
    // ============================================================================================

    TCanvas* c_rho_z_d_cosmo_before = new TCanvas("c_rho_z_d_cosmo_before", "Delayed vertex distribution (Cosmo before)", 1000, 1000);
    c_rho_z_d_cosmo_before->cd();

    h_rho_z_d_cosmo_before->Draw();

    c_rho_z_d_cosmo_before->Update();

    // ============================================================================================
    // Cosmo after - Prompt energy
    // ============================================================================================

    TCanvas* c_e_p_cosmo_after = new TCanvas("c_e_p_cosmo_after", "Prompt energy (Cosmo after)", 1000, 1000);
    c_e_p_cosmo_after->cd();

    h_e_p_cosmo_after->SetLineWidth(3);
    h_e_p_cosmo_after->SetLineStyle(kSolid);
    h_e_p_cosmo_after->SetLineColorAlpha(kBlue, 1.0);

    h_e_p_cosmo_after->Draw();

    c_e_p_cosmo_after->Update();

    // ============================================================================================
    // Cosmo after - Delayed energy
    // ============================================================================================

    TCanvas* c_e_d_cosmo_after = new TCanvas("c_e_d_cosmo_after", "Delayed energy (Cosmo after)", 1000, 1000);
    c_e_d_cosmo_after->cd();
    
    h_e_d_cosmo_after->SetLineWidth(3);
    h_e_d_cosmo_after->SetLineStyle(kSolid);
    h_e_d_cosmo_after->SetLineColorAlpha(kBlue, 1.0);

    h_e_d_cosmo_after->Draw();

    c_e_d_cosmo_after->Update();

    // ============================================================================================
    // Cosmo after - Prompt-Delayed time difference
    // ============================================================================================

    TCanvas* c_dt_cosmo_after = new TCanvas("c_dt_cosmo_after", "Prompt-Delayed time difference (Cosmo after)", 1000, 1000);
    c_dt_cosmo_after->cd();

    h_dt_cosmo_after->SetLineWidth(3);
    h_dt_cosmo_after->SetLineStyle(kSolid);
    h_dt_cosmo_after->SetLineColorAlpha(kBlue, 1.0);

    h_dt_cosmo_after->Draw();

    c_dt_cosmo_after->Update();

    // ============================================================================================
    // Cosmo after - Prompt-Delayed distance
    // ============================================================================================

    TCanvas* c_dr_cosmo_after = new TCanvas("c_dr_cosmo_after", "Prompt-Delayed distance (Cosmo after)", 1000, 1000);
    c_dr_cosmo_after->cd();

    h_dr_cosmo_after->Draw();

    c_dr_cosmo_after->Update();

    // ============================================================================================
    // Cosmo after - Prompt vertex position
    // ============================================================================================

    TCanvas* c_rho_z_p_cosmo_after = new TCanvas("c_rho_z_p_cosmo_after", "Prompt vertex distribution (Cosmo after)", 1000, 1000);
    c_rho_z_p_cosmo_after->cd();

    h_rho_z_p_cosmo_after->Draw();

    c_rho_z_p_cosmo_after->Update();

    // ============================================================================================
    // Cosmo after - Delayed vertex position
    // ============================================================================================

    TCanvas* c_rho_z_d_cosmo_after = new TCanvas("c_rho_z_d_cosmo_after", "Delayed vertex distribution (Cosmo after)", 1000, 1000);
    c_rho_z_d_cosmo_after->cd();

    h_rho_z_d_cosmo_after->Draw();

    c_rho_z_d_cosmo_after->Update();

    // ============================================================================================
    // Cosmo diff - Prompt energy
    // ============================================================================================

    TCanvas* c_e_p_cosmo_diff = new TCanvas("c_e_p_cosmo_diff", "Prompt energy (Cosmo diff)", 1000, 1000);
    c_e_p_cosmo_diff->cd();

    h_e_p_cosmo_diff->SetLineWidth(3);
    h_e_p_cosmo_diff->SetLineStyle(kSolid);
    h_e_p_cosmo_diff->SetLineColorAlpha(kBlue, 1.0);

    h_e_p_cosmo_diff->Draw();

    c_e_p_cosmo_diff->Update();

}