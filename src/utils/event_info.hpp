#ifndef ANALYSISGROUPC_UTILS_EVENTINFO_HPP_
#define ANALYSISGROUPC_UTILS_EVENTINFO_HPP_

#include <vector>

#include "event/IBD.hpp"

struct mult_info {
    vertex vtx;
    int type;
};

struct ibd_info {

    ibd pair;
    std::vector<vertex> neus;
    std::vector<mult_info> mults;

    ibd_info(const vertex& prompt, const vertex& delayed) :
        pair{prompt, delayed}
    {}

};

#endif // ANALYSISGROUPC_UTILS_EVENTINFO_HPP_