#ifndef UTILS_PLOT_HPP_
#define UTILS_PLOT_HPP_

#include <string>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>

#include "utils/event.hpp"
#include "utils/numpy.hpp"

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
    // double xmin = 0.7;
    // double xmax = 12.0;
    // double width = 0.20;
    // int nbins = std::round((xmax - xmin) / width) + 1;
    // std::vector<double> bins = np::linspace(xmin, xmax, nbins);
    std::vector<double> bins = create_custom_e_p_bins();
    TH1D* h = new TH1D(name.c_str(), title.c_str(), bins.size() - 1, bins.data());
    for (const ibd& i : ibds) {
        h->Fill(i.prompt.e);
    }
    return h;
}

TH1D* make_prompt_energy_plot(const std::string& name, const std::string& title, const std::vector<cosmogenic>& cosmos) {
    // double xmin = 0.7;
    // double xmax = 12.0;
    // double width = 0.20;
    // int nbins = std::round((xmax - xmin) / width) + 1;
    // std::vector<double> bins = np::linspace(xmin, xmax, nbins);
    std::vector<double> bins = create_custom_e_p_bins();
    TH1D* h = new TH1D(name.c_str(), title.c_str(), bins.size() - 1, bins.data());
    for (const cosmogenic& c : cosmos) {
        h->Fill(c.prompt.e);
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

void pimp_my_histogram(TH1D* h, Style_t linestyle, Width_t linewidth, Color_t linecolor, Float_t linealpha) {
    h->SetLineStyle(linestyle);
    h->SetLineWidth(linewidth);
    h->SetLineColorAlpha(linecolor, linealpha);
}

TCanvas* plot_basic(TH1D* h) {
    TCanvas* c = new TCanvas(Form("c_%s", h->GetName()), h->GetTitle(), 1000, 1000);
    c->cd();
    h->Draw();
    c->Update();
    return c;
}

#endif // UTILS_PLOT_HPP_