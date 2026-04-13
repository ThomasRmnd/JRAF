#ifndef UTILS_MUON_LOOKUP_HPP_
#define UTILS_MUON_LOOKUP_HPP_

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "utils/navigator.hpp"
#include "utils/timestamp.hpp"

class is_cd_muon_lookup {

public:

    void fill(const std::shared_ptr<basic_navigator>& nav) {
        m_times.clear();
        for (std::size_t k = 0ul; k < nav->method_mu.size(); ++k) {
            if (nav->totq_cd_mu[k] <= 0.0) continue;
            m_times.push_back(timestamp{nav->sec_mu[k], nav->nsec_mu[k]});
        }
        std::sort(m_times.begin(), m_times.end());
    }

    bool operator[](const timestamp& ts) const {
        if (m_times.empty()) return false;
        timestamp low_bound = ts - window;
        timestamp high_bound = ts + window;
        std::vector<timestamp>::const_iterator it_low = std::lower_bound(m_times.begin(), m_times.end(), low_bound);
        std::vector<timestamp>::const_iterator it_high = std::upper_bound(m_times.begin(), m_times.end(), high_bound);
        return std::distance(it_low, it_high) > 0l;
    }

private:

    std::vector<timestamp> m_times;
    const timestamp window{0, 1000};

};

class multiplicity_muon_lookup {

public:

    void fill(const std::shared_ptr<basic_navigator>& nav, const std::string& target) {
        m_times.clear();
        for (std::size_t k = 0ul; k < nav->method_mu.size(); ++k) {
            if (nav->method_mu[k] != target) continue;
            m_times.push_back(timestamp{nav->sec_mu[k], nav->nsec_mu[k]});
        }
        std::sort(m_times.begin(), m_times.end());
    }

    std::size_t operator[](const timestamp& ts) const {
        if (m_times.empty()) return 0ul;
        timestamp low_bound = ts - window;
        timestamp high_bound = ts + window;
        std::vector<timestamp>::const_iterator it_low = std::lower_bound(m_times.begin(), m_times.end(), low_bound);
        std::vector<timestamp>::const_iterator it_high = std::upper_bound(m_times.begin(), m_times.end(), high_bound);
        return std::distance(it_low, it_high);
    }

private:

    std::vector<timestamp> m_times;
    const timestamp window{0, 1000};

};

class stopping_muon_lookup {

public:

    void fill(const std::shared_ptr<basic_navigator>& nav, const std::string& target) {
        m_times.clear();
        const vec3 dummy_point{-999999.0, -999999.0, -999999.0};
        for (std::size_t k = 0ul; k < nav->method_mu.size(); ++k) {
            if (nav->method_mu[k] != target) continue;
            vec3 ipos{nav->posx_mu[k], nav->posy_mu[k], nav->posz_mu[k]};
            vec3 dir = unit(vec3{nav->dirx_mu[k], nav->diry_mu[k], nav->dirz_mu[k]});
            vec3 to_dummy = unit(dummy_point - ipos);
            double alpha = angle(dir, to_dummy);
            if (alpha > epsilon) continue;
            m_times.push_back(timestamp{nav->sec_mu[k], nav->nsec_mu[k]});
        }
        std::sort(m_times.begin(), m_times.end());
    }
        
    bool operator[](const timestamp& ts) const {
        if (m_times.empty()) return false;
        timestamp low_bound = ts - window;
        timestamp high_bound = ts + window;
        std::vector<timestamp>::const_iterator it_low = std::lower_bound(m_times.begin(), m_times.end(), low_bound);
        std::vector<timestamp>::const_iterator it_high = std::upper_bound(m_times.begin(), m_times.end(), high_bound);
        return std::distance(it_low, it_high) > 0l;
    }

private:

    std::vector<timestamp> m_times;
    const timestamp window{0, 1000};
    const double epsilon = 1e-4;

};

#endif // UTILS_MUON_LOOKUP_HPP_