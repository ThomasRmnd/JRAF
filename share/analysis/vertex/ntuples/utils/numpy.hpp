#ifndef UTILS_NUMPY_HPP_
#define UTILS_NUMPY_HPP_

#include <iostream>
#include <vector>

namespace np {

std::vector<double> linspace(double start, double stop, int num) {
    if (num <= 1) {
        if (num == 1) return {start};
        std::cerr << "Warning: linspace_cpp requires num >= 1. Returning empty vector\n";
        return {};
    }
    double step = (stop - start) / (num - 1);
    std::vector<double> result;
    for (int i = 0; i < num; ++i) {
        double value = start + static_cast<double>(i) * step;
        if (i == num - 1) {
            result.push_back(stop);
        } 
        else {
            result.push_back(value);
        }
    }
    
    return result;
}

} // namespace np

#endif // UTILS_NUMPY_HPP_