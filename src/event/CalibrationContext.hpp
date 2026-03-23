#ifndef JRAF_EVENT_CALIBRATIONCONTEXT_HPP_
#define JRAF_EVENT_CALIBRATIONCONTEXT_HPP_

#include <limits>

struct calibration_context {

    double totq = 0.0;
    double meanq = 0.0, stdq = 0.0;
    double minq = std::numeric_limits<double>::infinity(), maxq = 0.0;
    double meant = 0.0, stdt = 0.0;
    std::size_t npmt = 0ul;
    std::size_t nhit = 0ul;
    double meanhit = 0.0, stdhit = 0.0;

};

#endif // JRAF_EVENT_CALIBRATIONCONTEXT_HPP_