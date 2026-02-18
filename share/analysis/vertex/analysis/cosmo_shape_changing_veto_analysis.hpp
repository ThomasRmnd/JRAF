#ifndef ANALYSIS_COSMO_SHAPE_CHANGING_VETO_ANALYSIS_HPP_
#define ANALYSIS_COSMO_SHAPE_CHANGING_VETO_ANALYSIS_HPP_

#include <TLegend.h>

#include "analysis/cosmo_shape_analysis.hpp"
#include "utils/plot.hpp"

class cosmo_shape_changing_veto_analysis : public cosmo_shape_analysis {

public:

    using cosmo_shape_analysis::cosmo_shape_analysis;

    bool selection() override {
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

        if ( std::pow((m_nav->meta_prompt.stdhit - 0.55) / 0.45, 2.0) + std::pow((m_nav->meta_prompt.stdt - 170.0) / 80.0, 2.0) > 1.0 ) return false;

        m_dlat_mu2p.clear();
        m_dlat_mu2d.clear();
        m_dt_mu2p.clear();
        m_dt_mu2d.clear();
        m_is_sig.clear();

        multiplicity_muon_lookup nb_muons_in_cd_event;
        multiplicity_muon_lookup nb_muons_in_wp_event;
        nb_muons_in_cd_event.fill(m_nav, "CdClassify");
        nb_muons_in_wp_event.fill(m_nav, "WpBasic");

        // stopping_muon_lookup has_stopping_in_cd_event;
        stopping_muon_lookup has_stopping_in_wp_event;
        // has_stopping_in_cd_event.fill(m_nav, "CdClassify");
        has_stopping_in_wp_event.fill(m_nav, "WpBasic");

        for (std::size_t k = 0ul; k < m_nav->method_mu.size(); ++k) {
            if (m_nav->method_mu[k] != m_recname) continue;
            timestamp ts_mu{m_nav->sec_mu[k], m_nav->nsec_mu[k]};
            if (nb_muons_in_cd_event[ts_mu] > 1ul || nb_muons_in_wp_event[ts_mu] > 1ul) continue;
            if (has_stopping_in_wp_event[ts_mu]) continue;

            bool is_in_bkg = (
                ts_mu + m_ts_bkg_low <= m_nav->prompt.ts && m_nav->prompt.ts <= ts_mu + m_ts_bkg_high &&
                ts_mu + m_ts_bkg_low <= m_nav->delayed.ts && m_nav->delayed.ts <= ts_mu + m_ts_bkg_high
            );
            bool is_in_sig = (
                ts_mu + m_ts_sig_low <= m_nav->prompt.ts && m_nav->prompt.ts <= ts_mu + m_ts_sig_high &&
                ts_mu + m_ts_sig_low <= m_nav->delayed.ts && m_nav->delayed.ts <= ts_mu + m_ts_sig_high
            );
            vec3 pos_mu{m_nav->posx_mu[k], m_nav->posy_mu[k], m_nav->posz_mu[k]};
            vec3 dir_mu = unit(vec3{m_nav->dirx_mu[k], m_nav->diry_mu[k], m_nav->dirz_mu[k]});
            double d_mu2p = mag(cross(dir_mu, m_nav->prompt.pos - pos_mu));
            double d_mu2d = mag(cross(dir_mu, m_nav->delayed.pos - pos_mu));
            
            double radius_sig = 0.0;
            double radius_bkg = 0.0;
            if (is_in_sig) {
                radius_sig = m_radius + m_radius / timestamp_to_double(m_ts_sig_low - m_ts_sig_high) * timestamp_to_double(m_nav->prompt.ts - ts_mu);
            }
            if (is_in_bkg) {
                radius_bkg = m_radius + m_radius / timestamp_to_double(m_ts_bkg_high - m_ts_bkg_low) * timestamp_to_double(m_nav->prompt.ts - ts_mu);
            }

            if (radius_sig < d_mu2p && radius_bkg < d_mu2d) continue;
            m_dlat_mu2p.push_back(d_mu2p);
            m_dlat_mu2d.push_back(d_mu2d);
            m_dt_mu2p.push_back(timestamp_to_double(m_nav->prompt.ts - ts_mu));
            m_dt_mu2d.push_back(timestamp_to_double(m_nav->delayed.ts - ts_mu));
            m_is_sig.push_back(is_in_sig);
        }

        return (!m_is_sig.empty());
    }

};

#endif // ANALYSIS_COSMO_SHAPE_CHANGING_VETO_ANALYSIS_HPP_