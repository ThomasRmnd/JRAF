#ifndef ANALYSIS_COSMO_RATE_NEUTRON_VETO_ANALYSIS_HPP_
#define ANALYSIS_COSMO_RATE_NEUTRON_VETO_ANALYSIS_HPP_

#include "analysis/cosmo_rate_analysis.hpp"

class cosmo_rate_neutron_veto_analysis : public cosmo_rate_analysis {

    // Regarder nombre cosmo en fonction du nombre de neutron ==> est-ce que ça suit une loi de poisson

public:

    cosmo_rate_neutron_veto_analysis(const std::string& name, const std::string& filepath, const std::string& suffix) :
        cosmo_rate_analysis{name, filepath, suffix}
    {}

    virtual bool selection() override {
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

        return (nb_multu_veto == 0ul);
    }

};

#endif // ANALYSIS_COSMO_RATE_NEUTRON_VETO_ANALYSIS_HPP_