#ifndef UTILS_PLOT_HPP_
#define UTILS_PLOT_HPP_

#include <string>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TPaveStats.h>

#include "utils/event.hpp"
#include "utils/numpy.hpp"

enum class StatOpt {
    None        =         0,
    Name        =         1,  // n = 1
    Entries     =        10,  // e = 1
    Mean        =       100,  // m = 1
    MeanError   =       200,  // m = 2
    RMS         =      1000,  // r = 1
    RMSError    =      2000,  // r = 2
    Underflow   =     10000,  // u = 1
    Overflow    =    100000,  // o = 1
    Integral    =   1000000,  // i = 1
    Skew        =  10000000,  // s = 1
    SkewError   =  20000000,  // s = 2
    Kurtosis    = 100000000,  // k = 1
    KurtError   = 200000000,  // k = 2
};

inline StatOpt operator|(StatOpt a, StatOpt b) {
    return static_cast<StatOpt>(static_cast<int>(a) + static_cast<int>(b));
}

inline int to_root_opt(StatOpt mode) {
    return static_cast<int>(mode);
}

enum class FitOpt {
    None        =    0,
    Params      =    1,  // v = 1 (requires FitErrors)
    AllParams   =    2,  // v = 2
    Errors      =   10,  // e = 1
    Chi2NDF     =  100,  // c = 1
    Proba       = 1000,  // p = 1
};

inline FitOpt operator|(FitOpt a, FitOpt b) {
    return static_cast<FitOpt>(static_cast<int>(a) + static_cast<int>(b));
}

inline int to_root_opt(FitOpt mode) {
    return static_cast<int>(mode);
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

TH1D* make_prompt_energy_plot(const std::string& name, const std::string& title, const std::vector<ibd>& ibds) {
    std::vector<double> bins = create_custom_e_p_bins();
    TH1D* h = new TH1D(name.c_str(), title.c_str(), bins.size() - 1, bins.data());
    for (const ibd& i : ibds) {
        h->Fill(i.prompt.e);
    }
    return h;
}

TH1D* make_prompt_energy_plot(const std::string& name, const std::string& title, const std::vector<cosmogenic>& cosmos) {
    std::vector<double> bins = create_custom_e_p_bins();
    TH1D* h = new TH1D(name.c_str(), title.c_str(), bins.size() - 1, bins.data());
    for (const cosmogenic& c : cosmos) {
        h->Fill(c.prompt.e);
    }
    return h;
}

TH1D* make_normal_prompt_energy_plot(const std::string& name, const std::string& title, const std::vector<cosmogenic>& cosmos) {
    double xmin = 0.7;
    double xmax = 12.0;
    double width = 0.20;
    int nbins = std::round((xmax - xmin) / width) + 1;
    std::vector<double> bins = np::linspace(xmin, xmax, nbins);
    TH1D* h = new TH1D(name.c_str(), title.c_str(), bins.size() - 1, bins.data());
    for (const cosmogenic& c : cosmos) {
        h->Fill(c.prompt.e);
    }
    return h;
}

TH1D* make_normal_prompt_energy_plot(const std::string& name, const std::string& title, const std::vector<ibd>& ibds) {
    double xmin = 0.7;
    double xmax = 12.0;
    double width = 0.20;
    int nbins = std::round((xmax - xmin) / width) + 1;
    std::vector<double> bins = np::linspace(xmin, xmax, nbins);
    TH1D* h = new TH1D(name.c_str(), title.c_str(), bins.size() - 1, bins.data());
    for (const ibd& i : ibds) {
        h->Fill(i.prompt.e);
    }
    return h;
}

TH1D* make_delayed_energy_plot(const std::string& name, const std::string& title, const std::vector<cosmogenic>& cosmos) {
    double xmin = 2.0;
    double xmax = 2.5;
    double width = 0.02;
    int nbins = std::round((xmax - xmin) / width) + 1;
    std::vector<double> bins = np::linspace(xmin, xmax, nbins);
    TH1D* h = new TH1D(name.c_str(), title.c_str(), bins.size() - 1, bins.data());
    for (const cosmogenic& c : cosmos) {
        h->Fill(c.delayed.e);
    }
    return h;
}

TH1D* make_delayed_energy_plot(const std::string& name, const std::string& title, const std::vector<ibd>& ibds) {
    double xmin = 2.0;
    double xmax = 2.5;
    double width = 0.02;
    int nbins = std::round((xmax - xmin) / width) + 1;
    std::vector<double> bins = np::linspace(xmin, xmax, nbins);
    TH1D* h = new TH1D(name.c_str(), title.c_str(), bins.size() - 1, bins.data());
    for (const ibd& i : ibds) {
        h->Fill(i.delayed.e);
    }
    return h;
}

TH1D* make_prompt_delayed_time_plot(const std::string& name, const std::string& title, const std::vector<cosmogenic>& cosmos) {
    double xmin = 0.0;
    double xmax = 1.0;
    double width = 0.025;
    int nbins = std::round((xmax - xmin) / width) + 1;
    std::vector<double> bins = np::linspace(xmin, xmax, nbins);
    TH1D* h = new TH1D(name.c_str(), title.c_str(), bins.size() - 1, bins.data());
    for (const cosmogenic& c : cosmos) {
        h->Fill(timestamp_to_double(c.delayed.ts - c.prompt.ts) * 1000.0);
    }
    return h;
}

TH1D* make_prompt_delayed_time_plot(const std::string& name, const std::string& title, const std::vector<ibd>& ibds) {
    double xmin = 0.0;
    double xmax = 1.0;
    double width = 0.025;
    int nbins = std::round((xmax - xmin) / width) + 1;
    std::vector<double> bins = np::linspace(xmin, xmax, nbins);
    TH1D* h = new TH1D(name.c_str(), title.c_str(), bins.size() - 1, bins.data());
    for (const ibd& i : ibds) {
        h->Fill(timestamp_to_double(i.delayed.ts - i.prompt.ts) * 1000.0);
    }
    return h;
}

TH1D* make_prompt_delayed_distance_plot(const std::string& name, const std::string& title, const std::vector<cosmogenic>& cosmos) {
    double xmin = 0.0;
    double xmax = 1.5;
    double width = 0.05;
    int nbins = std::round((xmax - xmin) / width) + 1;
    std::vector<double> bins = np::linspace(xmin, xmax, nbins);
    TH1D* h = new TH1D(name.c_str(), title.c_str(), bins.size() - 1, bins.data());
    for (const cosmogenic& c : cosmos) {
        h->Fill(mag(c.delayed.pos - c.prompt.pos));
    }
    return h;
}

TH1D* make_prompt_delayed_distance_plot(const std::string& name, const std::string& title, const std::vector<ibd>& ibds) {
    double xmin = 0.0;
    double xmax = 1.5;
    double width = 0.05;
    int nbins = std::round((xmax - xmin) / width) + 1;
    std::vector<double> bins = np::linspace(xmin, xmax, nbins);
    TH1D* h = new TH1D(name.c_str(), title.c_str(), bins.size() - 1, bins.data());
    for (const ibd& i : ibds) {
        h->Fill(mag(i.delayed.pos - i.prompt.pos));
    }
    return h;
}

TH2D* make_prompt_spatial_plot(const std::string& name, const std::string& title, const std::vector<cosmogenic>& cosmos) {
    double xmin = 0.0;
    double xmax = 17.7 * 17.7;
    int xnbins = 51;
    double ymin = -20.0;
    double ymax = 20.0;
    int ynbins = 51;
    std::vector<double> xbins = np::linspace(xmin, xmax, xnbins);
    std::vector<double> ybins = np::linspace(ymin, ymax, ynbins);
    TH2D* h = new TH2D(name.c_str(), title.c_str(), xbins.size() - 1, xbins.data(), ybins.size() - 1, ybins.data());
    for (const cosmogenic& c : cosmos) {
        h->Fill((c.prompt.pos.x * c.prompt.pos.x + c.prompt.pos.y * c.prompt.pos.y) / 1.0e6, c.prompt.pos.z / 1000.0);
    }
    return h;
}

TH2D* make_prompt_spatial_plot(const std::string& name, const std::string& title, const std::vector<ibd>& ibds) {
    double xmin = 0.0;
    double xmax = 17.7 * 17.7;
    int xnbins = 51;
    double ymin = -20.0;
    double ymax = 20.0;
    int ynbins = 51;
    std::vector<double> xbins = np::linspace(xmin, xmax, xnbins);
    std::vector<double> ybins = np::linspace(ymin, ymax, ynbins);
    TH2D* h = new TH2D(name.c_str(), title.c_str(), xbins.size() - 1, xbins.data(), ybins.size() - 1, ybins.data());
    for (const ibd& i : ibds) {
        h->Fill((i.prompt.pos.x * i.prompt.pos.x + i.prompt.pos.y * i.prompt.pos.y) / 1.0e6, i.prompt.pos.z / 1000.0);
    }
    return h;
}

TH2D* make_delayed_spatial_plot(const std::string& name, const std::string& title, const std::vector<cosmogenic>& cosmos) {
    double xmin = 0.0;
    double xmax = 17.7 * 17.7;
    int xnbins = 51;
    double ymin = -20.0;
    double ymax = 20.0;
    int ynbins = 51;
    std::vector<double> xbins = np::linspace(xmin, xmax, xnbins);
    std::vector<double> ybins = np::linspace(ymin, ymax, ynbins);
    TH2D* h = new TH2D(name.c_str(), title.c_str(), xbins.size() - 1, xbins.data(), ybins.size() - 1, ybins.data());
    for (const cosmogenic& c : cosmos) {
        h->Fill((c.delayed.pos.x * c.delayed.pos.x + c.delayed.pos.y * c.delayed.pos.y) / 1.0e6, c.delayed.pos.z / 1000.0);
    }
    return h;
}

TH2D* make_delayed_spatial_plot(const std::string& name, const std::string& title, const std::vector<ibd>& ibds) {
    double xmin = 0.0;
    double xmax = 17.7 * 17.7;
    int xnbins = 51;
    double ymin = -20.0;
    double ymax = 20.0;
    int ynbins = 51;
    std::vector<double> xbins = np::linspace(xmin, xmax, xnbins);
    std::vector<double> ybins = np::linspace(ymin, ymax, ynbins);
    TH2D* h = new TH2D(name.c_str(), title.c_str(), xbins.size() - 1, xbins.data(), ybins.size() - 1, ybins.data());
    for (const ibd& i : ibds) {
        h->Fill((i.delayed.pos.x * i.delayed.pos.x + i.delayed.pos.y * i.delayed.pos.y) / 1.0e6, i.delayed.pos.z / 1000.0);
    }
    return h;
}

struct NameConfig {

    const char* name = "";
    const char* title = "";

};

void pimp_my_name(TNamed* n, const NameConfig& config) {
    n->SetName(config.name);
    n->SetTitle(config.title);
}

struct LineConfig {

    Style_t style = kSolid;
    Width_t width = 1.0;
    Color_t color = kBlack;
    Float_t alpha = 1.0;


};

void pimp_my_line(TAttLine* l, const LineConfig& config) {
    l->SetLineStyle(config.style);
    l->SetLineWidth(config.width);
    l->SetLineColorAlpha(config.color, config.alpha);
}

struct AxisConfig {

    struct Label {
        Style_t font = 42;
        Float_t size = 0.04;
        Float_t offset = 0.005;
        Color_t color = kBlack;
        Float_t alpha = 1.0;
    };

    struct Title {
        Style_t font = 42;
        Float_t size = 0.04;
        Float_t offset = 1.0;
        Color_t color = kBlack;
    };

    struct Tick {
        Float_t length = 0.03;
        Float_t size = 0.03;
    };

    Color_t color = kBlack;
    Int_t ndivisions = 510;
    Int_t maxdigits = 5;
    Label label;
    Title title;
    Tick tick;

};

void pimp_my_axis(TAttAxis* a, const AxisConfig& config) {
    a->SetAxisColor(config.color);
    a->SetNdivisions(config.ndivisions);
    a->SetMaxDigits(config.maxdigits);
    a->SetLabelFont(config.label.font);
    a->SetLabelSize(config.label.size);
    a->SetLabelOffset(config.label.offset);
    a->SetLabelColor(config.label.color);
    a->SetTitleFont(config.title.font);
    a->SetTitleSize(config.title.size);
    a->SetTitleOffset(config.title.offset);
    a->SetTitleColor(config.title.color);
    a->SetTickLength(config.tick.length);
    a->SetTickSize(config.tick.size);
}

struct PaveConfig {

};

TPaveStats* change_stats(TH1D* h, double xmin, double ymin, double xmax, double ymax, StatOpt statopt, FitOpt fitopt) {
    TPaveStats* st = (TPaveStats*)h->FindObject("stats");
    if (statopt == StatOpt::None && fitopt == FitOpt::None) {
        h->SetStats(false);
        return nullptr;
    }
    st->SetOptStat(to_root_opt(statopt));
    st->SetOptFit(to_root_opt(fitopt));
    st->SetX1NDC(xmin);
    st->SetX2NDC(xmax);
    st->SetY1NDC(ymin);
    st->SetY2NDC(ymax);
    return st;
}

TF1* create_exponential_decay_function(TH1D* h, double xmin, double xmax) {
    double constant_term = h->GetBinContent(h->GetNbinsX());
    double exponential_term = h->GetMaximum() - constant_term;
    double decay_term = h->GetRMS();
    TF1* f = new TF1(Form("f_%s", h->GetName()), "[0] * exp(-x / [1]) + [2]", xmin, xmax);
    f->SetParameter(0, exponential_term);
    f->SetParameter(1, decay_term);
    f->SetParameter(2, constant_term);
    return f;
}

TCanvas* plot_basic(TH1D* h, const char* options = "") {
    TCanvas* c = new TCanvas(Form("c_%s", h->GetName()), h->GetTitle(), 1000, 1000);
    c->cd();
    h->Draw(options);
    c->SetTickx();
    c->SetTicky();
    c->SetGrid();
    c->Update();
    return c;
}

TCanvas* plot_multiple(const std::string& name, const std::string& title, std::initializer_list<TH1D*> hists, TLegend* leg, const char* options = "") {
    if (hists.size() == 0) return nullptr;
    
    double max_val = 0;
    double min_val = std::numeric_limits<double>::max();
    for (TH1D* h : hists) {
        if (!h) continue;
        max_val = std::max(max_val, h->GetMaximum());
        min_val = std::min(min_val, h->GetMinimum());
    }
    
    for (TH1D* h : hists) {
        if (!h) continue;
        h->SetMaximum(max_val * 1.15);
        h->SetMinimum(std::min(0.0, min_val));
    }
    
    TCanvas* c = new TCanvas(name.c_str(), title.c_str(), 1000, 1000);
    c->SetGrid();
    c->cd();
    
    leg->SetBorderSize(1);
    
    bool is_first = true;
    for (TH1D* h : hists) {
        if (!h) continue;
        if (is_first) {
            h->Draw(options);
            is_first = false;
        } 
        else {
            h->Draw(Form("%s SAME", options));
        }
        leg->AddEntry(h, h->GetTitle(), "l");
    }
    
    leg->Draw();

    c->SetTickx();
    c->SetTicky();
    c->SetGrid();
    c->Update();
    return c;
}

TCanvas* plot_basic(TH2D* h, const char* options = "") {
    TCanvas* c = new TCanvas(Form("c_%s", h->GetName()), h->GetTitle(), 1000, 1000);
    c->cd();
    h->Draw(options);
    c->SetTickx();
    c->SetTicky();
    c->SetGrid();
    c->Update();
    return c;
}

#endif // UTILS_PLOT_HPP_