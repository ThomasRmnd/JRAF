#ifndef ANLYSIS_IBD_MUON_VETO_ANALYSIS_HPP_
#define ANLYSIS_IBD_MUON_VETO_ANALYSIS_HPP_

#include "analysis/basic_analysis.hpp"

class ibd_muon_veto_analysis : public basic_analysis {

public:

    ibd_muon_veto_analysis(
        const std::string& name, 
        const std::string& filepath, const std::string& suffix,
        const std::string& recname,
        const timestamp& low, const timestamp& high, double radius
    ) :
        basic_analysis{name, filepath, suffix},
        m_recname{recname},
        m_ts_low{low},
        m_ts_high{high},
        m_radius{radius}
    {}

    bool selection() override {
        if (m_nav->meta_prompt.stdt > 200.0 || m_nav->meta_delayed.stdt > 200.0) return false; // Flasher cut
        if (mag(m_nav->prompt.pos) > 16500.0) return false; // Fiducial cut
        if ((m_nav->prompt.pos.z < -15500.0 || 15500 < m_nav->prompt.pos.z) && std::sqrt(m_nav->prompt.pos.x * m_nav->prompt.pos.x + m_nav->prompt.pos.y * m_nav->prompt.pos.y) < 3000.0) return false; // Chimney cut
        if (m_nav->prompt.e < 0.7 || 12.0 < m_nav->prompt.e) return false; // Prompt energy cut
        if (m_nav->delayed.e < 2.0 || 2.5 < m_nav->delayed.e) return false; // Delayed energy cut

        std::size_t nb_multu_veto = 0ul;
        for (std::size_t k = 0ul; k < m_nav->e_mult.size(); ++k) {
            if (m_nav->e_mult[k] < 2.0 || 12.0 < m_nav->e_mult[k]) continue;
            timestamp ts_mult{m_nav->sec_mult[k], m_nav->nsec_mult[k]};
            vec3 pos_mult{m_nav->posx_mult[k], m_nav->posy_mult[k], m_nav->posz_mult[k]};
            // if (pos_mult.Mag() > 17700.0) continue;
            // if ((pos_p - pos_n).Mag() > 4000.0 || (pos_p - pos_n).Mag() > 4000.0) continue;
            if (ts_mult < m_nav->prompt.ts - timestamp{0, 1000000} || m_nav->delayed.ts + timestamp{0, 1000000} < ts_mult) continue;
            ++nb_multu_veto;
        }
        if (nb_multu_veto != 0ul) return false;

        std::size_t nb_neutron_veto = 0ul;
        for (std::size_t k = 0ul; k < m_nav->e_n.size(); ++k) {
            if (m_nav->e_n[k] < 1.5 || 20.0 < m_nav->e_n[k]) continue;
            timestamp ts_n{m_nav->sec_n[k], m_nav->nsec_n[k]};
            vec3 pos_n{m_nav->posx_n[k], m_nav->posy_n[k], m_nav->posz_n[k]};
            // if (pos_n.Mag() > 17700.0) continue;
            if (mag(m_nav->prompt.pos - pos_n) > 4000.0 || mag(m_nav->delayed.pos - pos_n) > 4000.0) continue;
            if (m_nav->prompt.ts < ts_n + timestamp{0, 20000} || ts_n + timestamp{0, 1200000000} < m_nav->prompt.ts) continue;
            if (m_nav->delayed.ts < ts_n + timestamp{0, 20000} || ts_n + timestamp{0, 1200000000} < m_nav->delayed.ts) continue;
            ++nb_neutron_veto;
        }
        if (nb_neutron_veto != 0ul) return false;

        std::size_t nb_muon_veto = 0ul;
        for (std::size_t k = 0ul; k < m_nav->method_mu.size(); ++k) {
            if (m_nav->method_mu[k] != m_recname) continue;
            timestamp ts_mu{m_nav->sec_mu[k], m_nav->nsec_mu[k]};
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
        if (nb_muon_veto != 0ul) return false;

        return true;
    }

private:

    std::string m_recname;
    timestamp m_ts_low;
    timestamp m_ts_high;
    double m_radius;

};  

#endif // ANLYSIS_IBD_MUON_VETO_ANALYSIS_HPP_