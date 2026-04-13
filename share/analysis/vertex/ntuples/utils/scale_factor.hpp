#ifndef UTILS_SCALE_FACTOR_HPP_
#define UTILS_SCALE_FACTOR_HPP_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "utils/timestamp.hpp"

class TimeCorrector {

public:

    struct Point {
        timestamp ts;
        double factor;
    };

    typedef std::vector<Point>                              container_type;
    typedef typename container_type::size_type              size_type;
    typedef typename container_type::difference_type        difference_type;
    typedef typename container_type::const_reference        const_reference;
    typedef typename container_type::const_iterator         const_iterator;
    typedef typename container_type::const_reverse_iterator const_reverse_iterator;

    bool load(const std::string& filepath) {
        std::ifstream ifs(filepath);
        if (!ifs.is_open()) {
            std::cerr << "Cannot open file " << filepath << " for the time correction\n";
            return false;
        }
        std::string line, ts_str, factor_str;
        std::getline(ifs, line);
        while (std::getline(ifs, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            if (std::getline(ss, ts_str, ',') && std::getline(ss, factor_str, ',')) {
                long sec = std::stoll(ts_str);
                timestamp ts{static_cast<time_t>(sec), 0};
                data.push_back({ts, std::stod(factor_str)});
            }
        }
        std::sort(data.begin(), data.end());
        return true;
    }

    size_type size() const { return data.size(); }

    const_reference front() const { return data.front(); }
    const_reference back() const { return data.back(); }

    const_iterator begin() const { return data.begin(); }
    const_iterator cbegin() const { return data.cbegin(); }
    const_iterator end() const { return data.end(); }
    const_iterator cend() const { return data.cend(); }

    const_reverse_iterator rbegin() const { return data.rbegin(); }
    const_reverse_iterator crbegin() const { return data.crbegin(); }
    const_reverse_iterator rend() const { return data.rend(); }
    const_reverse_iterator crend() const { return data.crend(); }

    double interpolate(const timestamp& ts) const {
        if (data.empty()) return 1.0;
        std::vector<Point>::const_iterator it = std::lower_bound(data.begin(), data.end(), Point{ts, 0.0});
        if (it == data.begin()) return data.front().factor;
        if (it == data.end()) return data.back().factor;
        
        const Point& p1 = *(it - 1);
        const Point& p2 = *it;

        double fraction = timestamp_to_double(ts - p1.ts) / timestamp_to_double(p2.ts - p1.ts);
        return p1.factor + fraction * (p2.factor - p1.factor);
    }

private:

    std::vector<Point> data;

};

inline bool operator<(const TimeCorrector::Point& lhs, const TimeCorrector::Point& rhs) {
    return lhs.ts < rhs.ts;
}

class GlobalTimeCorrector {

public:

    bool load() {
        if (!m_tc_p25c.load("/cvmfs/juno.ihep.ac.cn/dbdata/main/dbdata/offline-data/Reconstruction/OMILREC/RecMap/nPEMap/Final_time_correction_P25C_AvgSPN.csv")) return false;
        if (!m_tc_p25d.load("/cvmfs/juno.ihep.ac.cn/dbdata/main/dbdata/offline-data/Reconstruction/OMILREC/RecMap/nPEMap/Final_time_correction_P25D_AvgSPN.csv")) return false;
        if (!m_tc_p25c.size() || !m_tc_p25d.size()) {
            std::cerr << "P25C or P25D time correction is empty\n";
            return false;
        }
        return true;
    }

    double interpolate(const timestamp& ts) {
        if (!m_tc_p25c.size() || !m_tc_p25d.size()) return 1.0;
        if (ts <= m_tc_p25c.back().ts) return m_tc_p25c.interpolate(ts);
        if (ts <= m_tc_p25d.back().ts) return m_tc_p25d.interpolate(ts);
        return 1.0;
    }

private:

    const int c_lower_run_id_p25c = 9789;
    const int c_upper_run_id_p25c = 11039;
    const int c_lower_run_id_p25d = 11049;
    const int c_upper_run_id_p25d = 12135;

    TimeCorrector m_tc_p25c;
    TimeCorrector m_tc_p25d;

};

#endif // UTILS_SCALE_FACTOR_HPP_