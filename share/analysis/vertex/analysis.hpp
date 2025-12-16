#ifndef UTILS_ANALYSIS_HPP_
#define UTILS_ANALYSIS_HPP_

#include <iostream>
#include <string>

#include <TChain.h>

#include "utils/event.hpp"
#include "utils/timestamp.hpp"

class analysis_base {

public:

    analysis_base(const std::string& name_) : name{name_} {}

    virtual ~analysis_base() = default;

    virtual bool retrieve(const std::string& filename) {
        m_chain = new TChain(name.c_str());
        if (!m_chain) {
            std::cerr << "Cannot create TChain " << name << '\n';
            return false;
        }
        m_chain->Add(filename.c_str());
        if (m_chain->GetEntries() == 0) {
            std::cerr << "No entries in TChain " << name << " from file " << filename << '\n';
            delete m_chain;
            m_chain = nullptr;
            return false;
        }

        m_chain->SetBranchAddress("run_id", &run_id);

        m_chain->SetBranchAddress("posx_p", &prompt.pos.x);
        m_chain->SetBranchAddress("posy_p", &prompt.pos.y);
        m_chain->SetBranchAddress("posz_p", &prompt.pos.z);
        m_chain->SetBranchAddress("e_p", &prompt.e);
        m_chain->SetBranchAddress("totq_p", &prompt.q);
        m_chain->SetBranchAddress("sec_p", &prompt.ts.sec);
        m_chain->SetBranchAddress("nsec_p", &prompt.ts.nsec);

        m_chain->SetBranchAddress("meanq_p", &meta_prompt.meanq);
        m_chain->SetBranchAddress("stdq_p", &meta_prompt.stdq);
        m_chain->SetBranchAddress("minq_p", &meta_prompt.minq);
        m_chain->SetBranchAddress("maxq_p", &meta_prompt.maxq);
        m_chain->SetBranchAddress("nhit_p", &meta_prompt.nhit);
        m_chain->SetBranchAddress("meant_p", &meta_prompt.meant);
        m_chain->SetBranchAddress("stdt_p", &meta_prompt.stdt);

        m_chain->SetBranchAddress("posx_d", &delayed.pos.x);
        m_chain->SetBranchAddress("posy_d", &delayed.pos.y);
        m_chain->SetBranchAddress("posz_d", &delayed.pos.z);
        m_chain->SetBranchAddress("e_d", &delayed.e);
        m_chain->SetBranchAddress("totq_d", &delayed.q);
        m_chain->SetBranchAddress("sec_d", &delayed.ts.sec);
        m_chain->SetBranchAddress("nsec_d", &delayed.ts.nsec);

        m_chain->SetBranchAddress("meanq_d", &meta_delayed.meanq);
        m_chain->SetBranchAddress("stdq_d", &meta_delayed.stdq);
        m_chain->SetBranchAddress("minq_d", &meta_delayed.minq);
        m_chain->SetBranchAddress("maxq_d", &meta_delayed.maxq);
        m_chain->SetBranchAddress("nhit_d", &meta_delayed.nhit);
        m_chain->SetBranchAddress("meant_d", &meta_delayed.meant);
        m_chain->SetBranchAddress("stdt_d", &meta_delayed.stdt);

        return true;
    }

    std::size_t size() const {
        if (!m_chain) return 0ul;
        return static_cast<std::size_t>(m_chain->GetEntries());
    }

    bool entry(std::ptrdiff_t n) {
        if (!m_chain) return false;
        if (n < 0 || static_cast<std::size_t>(n) >= size()) return false;
        m_chain->GetEntry(static_cast<Long64_t>(n));
        return true;
    }

    virtual bool selection() = 0;

    virtual void print() const {
        std::cout << "Prompt: " << prompt << '\n';
        std::cout << "Delayed: " << delayed << '\n';
    }

    std::string name;

    int run_id;

    vertex prompt;
    vertex delayed;

    vertex_metadata meta_prompt;
    vertex_metadata meta_delayed;

protected:

    TChain* m_chain = nullptr;

};

class basic_analysis : public analysis_base {

public:

    basic_analysis(const std::string& suffix) : analysis_base{"IBDAnalysis" + suffix} {}

    ~basic_analysis() override = default;

    bool retrieve(const std::string& filename) override {
        if (!analysis_base::retrieve(filename)) return false;
        
        m_chain->SetBranchAddress("posx_n", &posx_n);
        m_chain->SetBranchAddress("posy_n", &posy_n);
        m_chain->SetBranchAddress("posz_n", &posz_n);
        m_chain->SetBranchAddress("e_n", &e_n);
        m_chain->SetBranchAddress("totq_n", &totq_n);
        m_chain->SetBranchAddress("sec_n", &sec_n);
        m_chain->SetBranchAddress("nsec_n", &nsec_n);

        m_chain->SetBranchAddress("posx_mult", &posx_mult);
        m_chain->SetBranchAddress("posy_mult", &posy_mult);
        m_chain->SetBranchAddress("posz_mult", &posz_mult);
        m_chain->SetBranchAddress("e_mult", &e_mult);
        m_chain->SetBranchAddress("totq_mult", &totq_mult);
        m_chain->SetBranchAddress("sec_mult", &sec_mult);
        m_chain->SetBranchAddress("nsec_mult", &nsec_mult);
        m_chain->SetBranchAddress("mult_type", &mult_type);

        m_chain->SetBranchAddress("method_mu", &method_mu);
        m_chain->SetBranchAddress("loc_mu", &loc_mu);
        m_chain->SetBranchAddress("posx_mu", &posx_mu);
        m_chain->SetBranchAddress("posy_mu", &posy_mu);
        m_chain->SetBranchAddress("posz_mu", &posz_mu);
        m_chain->SetBranchAddress("dirx_mu", &dirx_mu);
        m_chain->SetBranchAddress("diry_mu", &diry_mu);
        m_chain->SetBranchAddress("dirz_mu", &dirz_mu);
        m_chain->SetBranchAddress("totq_mu", &totq_mu);
        m_chain->SetBranchAddress("sec_mu", &sec_mu);
        m_chain->SetBranchAddress("nsec_mu", &nsec_mu);
        m_chain->SetBranchAddress("quality_mu", &quality_mu);

        return true;
    }

    virtual bool selection() override {
        if (meta_prompt.stdt > 200.0 || meta_delayed.stdt > 200.0) return false; // Flasher cut
        if (mag(prompt.pos) > 16500.0) return false; // Fiducial cut
        if ((prompt.pos.z < -15500.0 || 15500 < prompt.pos.z) && std::sqrt(prompt.pos.x * prompt.pos.x + prompt.pos.y * prompt.pos.y) < 3000.0) return false; // Chimney cut
        if (prompt.e < 0.7 || 12.0 < prompt.e) return false; // Prompt energy cut
        if (delayed.e < 2.0 || 2.5 < delayed.e) return false; // Delayed energy cut

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < e_mult->size(); ++k) {
            if (e_mult->operator[](k) < 2.0 || 12.0 < e_mult->operator[](k)) continue;
            timestamp ts_mult{sec_mult->operator[](k), nsec_mult->operator[](k)};
            vec3 pos_mult{posx_mult->operator[](k), posy_mult->operator[](k), posz_mult->operator[](k)};
            // if (pos_mult.Mag() > 17700.0) continue;
            // if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_mult < prompt.ts - timestamp{0, 1000000} || delayed.ts + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }

        std::size_t nb_neutron_veto = 0ul;
        for (std::size_t k = 0ul; k < e_n->size(); ++k) {
            if (e_n->operator[](k) < 1.5 || 20.0 < e_n->operator[](k)) continue;
            timestamp ts_n{sec_n->operator[](k), nsec_n->operator[](k)};
            vec3 pos_n{posx_n->operator[](k), posy_n->operator[](k), posz_n->operator[](k)};
            // if (pos_n.Mag() > 17700.0) continue;
            if (mag(prompt.pos - pos_n) > 4000.0 || mag(delayed.pos - pos_n) > 4000.0) continue;
            if (prompt.ts < ts_n + timestamp{0, 20000} || ts_n + timestamp{0, 1200000000} < prompt.ts) continue;
            if (delayed.ts < ts_n + timestamp{0, 20000} || ts_n + timestamp{0, 1200000000} < delayed.ts) continue;
            ++nb_neutron_veto;
        }

        return (nb_neutron_veto == 0ul && nb_multu_veto == 0ul);
    }

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

};

using ibd_analysis = basic_analysis;

class cosmo_shape_analysis : public basic_analysis {

public:

    cosmo_shape_analysis(const std::string& suffix, const timestamp& low, const timestamp& high, double radius) :
        basic_analysis{suffix},
        m_ts_low{low},
        m_ts_high{high},
        m_radius{radius}
    {}

    bool selection() override {
        if (meta_prompt.stdt > 200.0 || meta_delayed.stdt > 200.0) return false; // Flasher cut
        if (mag(prompt.pos) > 16500.0) return false; // Fiducial cut
        if ((prompt.pos.z < -15500.0 || 15500 < prompt.pos.z) && std::sqrt(prompt.pos.x * prompt.pos.x + prompt.pos.y * prompt.pos.y) < 3000.0) return false; // Chimney cut
        if (prompt.e < 0.7 || 12.0 < prompt.e) return false; // Prompt energy cut
        if (delayed.e < 2.0 || 2.5 < delayed.e) return false; // Delayed energy cut

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < e_mult->size(); ++k) {
            if (e_mult->operator[](k) < 2.0 || 12.0 < e_mult->operator[](k)) continue;
            timestamp ts_mult{sec_mult->operator[](k), nsec_mult->operator[](k)};
            vec3 pos_mult{posx_mult->operator[](k), posy_mult->operator[](k), posz_mult->operator[](k)};
            // if (pos_mult.Mag() > 17700.0) continue;
            // if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_mult < prompt.ts - timestamp{0, 1000000} || delayed.ts + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }

        bool found = false;
        for (std::size_t k = 0ul; k < method_mu->size() && !found; ++k) {
            if (method_mu->operator[](k) != "CdWpTtChi2") continue;
            timestamp ts_mu{sec_mu->operator[](k), nsec_mu->operator[](k)};
            if (prompt.ts < ts_mu + m_ts_low || m_ts_high + ts_mu < prompt.ts) continue;
            if (delayed.ts < ts_mu + m_ts_low || m_ts_high + ts_mu < delayed.ts) continue;
            vec3 pos_mu{posx_mu->operator[](k), posy_mu->operator[](k), posz_mu->operator[](k)};
            vec3 dir_mu{dirx_mu->operator[](k), diry_mu->operator[](k), dirz_mu->operator[](k)};
            double d_mu2p = mag(cross(dir_mu, prompt.pos - pos_mu));
            double d_mu2d = mag(cross(dir_mu, delayed.pos - pos_mu));
            if (m_radius < d_mu2p || m_radius < d_mu2d) continue;
            found = true;
            dlat_mu2p = d_mu2p;
            dlat_mu2d = d_mu2d;
            dt_mu2p = timestamp_to_double(prompt.ts - ts_mu);
            dt_mu2d = timestamp_to_double(delayed.ts - ts_mu);
        }

        return (nb_multu_veto == 0ul && found);
    }

    double dlat_mu2p;
    double dlat_mu2d;
    double dt_mu2p;
    double dt_mu2d;

private:

    timestamp m_ts_low;
    timestamp m_ts_high;
    double m_radius;

};

class cosmo_rate_analysis : public basic_analysis {

public:

    cosmo_rate_analysis(const std::string& suffix) :
        basic_analysis{suffix}
    {}

    bool selection() override {
        if (meta_prompt.stdt > 200.0 || meta_delayed.stdt > 200.0) return false; // Flasher cut
        if (mag(prompt.pos) > 16500.0) return false; // Fiducial cut
        if ((prompt.pos.z < -15500.0 || 15500 < prompt.pos.z) && std::sqrt(prompt.pos.x * prompt.pos.x + prompt.pos.y * prompt.pos.y) < 3000.0) return false; // Chimney cut
        if (prompt.e < 0.7 || 12.0 < prompt.e) return false; // Prompt energy cut
        if (delayed.e < 2.0 || 2.5 < delayed.e) return false; // Delayed energy cut

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < e_mult->size(); ++k) {
            if (e_mult->operator[](k) < 2.0 || 12.0 < e_mult->operator[](k)) continue;
            timestamp ts_mult{sec_mult->operator[](k), nsec_mult->operator[](k)};
            vec3 pos_mult{posx_mult->operator[](k), posy_mult->operator[](k), posz_mult->operator[](k)};
            // if (pos_mult.Mag() > 17700.0) continue;
            // if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_mult < prompt.ts - timestamp{0, 1000000} || delayed.ts + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }

        return (nb_multu_veto == 0ul);
    }

};

#endif // UTILS_ANALYSIS_HPP_