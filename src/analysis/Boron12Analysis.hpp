#ifndef JRAF_ANALYSIS_BORON12ANALYSIS_HPP_
#define JRAF_ANALYSIS_BORON12ANALYSIS_HPP_

#include "analysis/Analysis.hpp"

class Boron12Analysis : public Analysis {

public:

    Boron12Analysis(const std::string& name, const std::string& method);

    ~Boron12Analysis() override = default;

    bool initialize() override;
    void process(const EventContext::View& events) override;

private:

    int run_id;

    double posx;
    double posy;
    double posz;
    double e;
    time_t sec;
    int nsec;

    double totq;
    double meanq;
    double stdq;
    double minq;
    double maxq;
    double meant;
    double stdt;
    std::size_t npmt;
    std::size_t nhit;
    double meanhit;
    double stdhit;

};

#endif // JRAF_ANALYSIS_BORON12ANALYSIS_HPP_