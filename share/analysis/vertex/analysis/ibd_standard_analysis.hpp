#ifndef ANALYSIS_IBD_STANDARD_ANALYSIS_HPP_
#define ANALYSIS_IBD_STANDARD_ANALYSIS_HPP_

#include <set>

#include <TFile.h>

#include "analysis/ibd_analysis.hpp"

class ibd_standard_analysis : public ibd_analysis {

public:

    using ibd_analysis::ibd_analysis;

    virtual ~ibd_standard_analysis() override = default;

    virtual bool selection() override {
        double e_p = m_nav->prompt.e / m_gtc.interpolate(m_nav->prompt.ts);
        double e_d = m_nav->delayed.e / m_gtc.interpolate(m_nav->delayed.ts);

        if (e_p < 0.7 || 12.0 < e_p) return false;
        if (e_d < 2.0 || 2.5 < e_d) return false;
        if (mag(m_nav->prompt.pos) > 16500.0) return false;
        if (std::abs(m_nav->prompt.pos.z) > 15500.0 && std::sqrt(m_nav->prompt.pos.x * m_nav->prompt.pos.x + m_nav->prompt.pos.y * m_nav->prompt.pos.y) < 2000.0) return false;
        timestamp ts_diff = m_nav->delayed.ts - m_nav->prompt.ts;
        if (ts_diff < timestamp{0, 5000} || timestamp{0, 1000000} < ts_diff) return false;
        vec3 pos_diff = m_nav->delayed.pos - m_nav->prompt.pos;
        if (mag(pos_diff) > 1500.0) return false;

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < m_nav->e_mult.size(); ++k) {
            timestamp ts_mult{m_nav->sec_mult[k], m_nav->nsec_mult[k]};
            vec3 pos_mult{m_nav->posx_mult[k], m_nav->posy_mult[k], m_nav->posz_mult[k]};
            double e_mult = m_nav->e_mult[k] / m_gtc.interpolate(ts_mult);
            if (e_mult < 2.0 || 12.0 < e_mult) continue;
            if (ts_mult < m_nav->prompt.ts - timestamp{0, 1000000} || m_nav->delayed.ts + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }
        if (nb_multu_veto) return false;

        std::size_t nb_neutron_veto = 0ul;
        for (std::size_t k = 0ul; k < m_nav->e_n.size(); ++k) {
            timestamp ts_n{m_nav->sec_n[k], m_nav->nsec_n[k]};
            vec3 pos_n{m_nav->posx_n[k], m_nav->posy_n[k], m_nav->posz_n[k]};
            double e_n = m_nav->e_n[k] / m_gtc.interpolate(ts_n);
            if (e_n < 1.5 || 20.0 < e_n) continue;
            if (m_nav->stdt_n[k] > 275.0) continue;
            if (
                (mag(m_nav->prompt.pos - pos_n) < 4000.0 && ts_n + timestamp{0, 20000} < m_nav->prompt.ts && m_nav->prompt.ts < ts_n + timestamp{0, 1200000000}) ||
                (mag(m_nav->delayed.pos - pos_n) < 4000.0 && ts_n + timestamp{0, 20000} < m_nav->delayed.ts && m_nav->delayed.ts < ts_n + timestamp{0, 1200000000})
            ) {
                ++nb_neutron_veto;
            }
        }
        if (!nb_neutron_veto) return false;

        if ( std::pow((m_nav->meta_prompt.stdhit - 0.55) / 0.45, 2.0) + std::pow((m_nav->meta_prompt.stdt - 170.0) / 80.0, 2.0) > 1.0 ) return false;

        calculate_dt_to_last_muon();
        calculate_dlat_dt_muon_to_prompt();

        return true;
    }

};

#endif // ANALYSIS_IBD_STANDARD_ANALYSIS_HPP_