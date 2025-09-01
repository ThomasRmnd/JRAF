#ifndef ANALYSISGROUPC_EVENT_IBD_HPP_
#define ANALYSISGROUPC_EVENT_IBD_HPP_

#include "event/Vertex.hpp"

struct ibd {

    vertex prompt;
    vertex delayed;

    ibd(const vertex& prompt_, const vertex& delayed_);

};

template<class _Char, class _Traits>
std::basic_ostream<_Char, _Traits>& operator<<(std::basic_ostream<_Char, _Traits>& os, const ibd& ibd) {
    return os << "prompt: {" << ibd.prompt << "}, delayed: {" << ibd.delayed << "}";
}

#endif // ANALYSISGROUPC_EVENT_IBD_HPP_