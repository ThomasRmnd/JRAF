#ifndef ANALYSIS_IBD_NO_NEUTRON_VETO_MUON_VETO_ANALYSIS_HPP_
#define ANALYSIS_IBD_NO_NEUTRON_VETO_MUON_VETO_ANALYSIS_HPP_

#include "analysis/ibd_analysis.hpp"

class ibd_no_neutron_veto_muon_veto_analysis : public ibd_analysis {

public:

    ibd_no_neutron_veto_muon_veto_analysis(
        const std::string& name, 
        const std::string& filepath, const std::string& suffix,
        const std::string& recname,
        const timestamp& low, const timestamp& high, double radius
    ) :
        ibd_analysis{name, filepath, suffix},
        m_recname{recname},
        m_ts_low{low},
        m_ts_high{high},
        m_radius{radius}
    {}

    ~ibd_no_neutron_veto_muon_veto_analysis() override = default;

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

        multiplicity_muon_lookup nb_muons_in_cd_event;
        multiplicity_muon_lookup nb_muons_in_wp_event;
        nb_muons_in_cd_event.fill(m_nav, "CdClassify");
        nb_muons_in_wp_event.fill(m_nav, "WpBasic");

        // stopping_muon_lookup has_stopping_in_cd_event;
        stopping_muon_lookup has_stopping_in_wp_event;
        // has_stopping_in_cd_event.fill(m_nav, "CdClassify");
        has_stopping_in_wp_event.fill(m_nav, "WpBasic");

        std::size_t nb_muon_veto = 0ul;
        for (std::size_t k = 0ul; k < m_nav->method_mu.size(); ++k) {
            if (m_nav->method_mu[k] != m_recname) continue;
            timestamp ts_mu{m_nav->sec_mu[k], m_nav->nsec_mu[k]};
            if (nb_muons_in_cd_event[ts_mu] > 1ul || nb_muons_in_wp_event[ts_mu] > 1ul) continue;
            if (has_stopping_in_wp_event[ts_mu]) continue;

            bool found_neutron = false;
            for (std::size_t l = 0ul; l < m_nav->e_n.size() && !found_neutron; ++l) {
                timestamp ts_n{m_nav->sec_n[l], m_nav->nsec_n[l]};
                if (ts_n < ts_mu + timestamp{0, 20000} || ts_mu + timestamp{0, 2000000} < ts_n) continue;
                found_neutron = true;
            }
            if (!found_neutron) continue;

            vec3 pos_mu{m_nav->posx_mu[k], m_nav->posy_mu[k], m_nav->posz_mu[k]};
            vec3 dir_mu{m_nav->dirx_mu[k], m_nav->diry_mu[k], m_nav->dirz_mu[k]};
            
            bool is_in_ts_veto = (
                ts_mu + m_ts_low < m_nav->prompt.ts && m_nav->prompt.ts < ts_mu + m_ts_high &&
                ts_mu + m_ts_low < m_nav->delayed.ts && m_nav->delayed.ts < ts_mu + m_ts_high
            );
            bool is_in_pos_veto = (
                mag(cross(dir_mu, m_nav->prompt.pos - pos_mu)) < m_radius &&
                mag(cross(dir_mu, m_nav->delayed.pos - pos_mu)) < m_radius
            );
            
            if (!is_in_ts_veto || !is_in_pos_veto) continue;
            ++nb_muon_veto;
        }
        if (nb_muon_veto) return false;

        if ( std::pow((m_nav->meta_prompt.stdhit - 0.55) / 0.45, 2.0) + std::pow((m_nav->meta_prompt.stdt - 170.0) / 80.0, 2.0) > 1.0 ) return false;

        calculate_dt_to_last_muon();
        calculate_dlat_dt_muon_to_prompt();

        return true;
    }

private:

    std::string m_recname;
    timestamp m_ts_low;
    timestamp m_ts_high;
    double m_radius;

};

#endif // ANALYSIS_IBD_NO_NEUTRON_VETO_MUON_VETO_ANALYSIS_HPP_