#ifndef JRAF_MUON_MUONTAGGER_HPP_
#define JRAF_MUON_MUONTAGGER_HPP_

#include "SniperKernel/ToolBase.h"

#include "Context/TimeStamp.h"
#include "EvtNavigator/NavBuffer.h"

#include "muon/TaggingContext.hpp"
#include "muon/TaggingResult.hpp"

class MuonTagger : public ToolBase {

public:

    MuonTagger(const std::string& name);

    virtual ~MuonTagger() override = default;

    virtual bool initialize() override;
    virtual muon_tagging_result tag(const muon_tagging_context& ctx);

protected:

    double m_cd_muon_totq_thold = 1000.0;
    double m_cd_only_muon_totq_thold = 30000.0;
    double m_wp_muon_totq_thold = 400.0;
    double m_wp_only_muon_totq_thold = 700.0;
    long long m_cd_afterpulse_thold_ns = 50000;
    long long m_wp_afterpulse_thold_ns = 4000;
    TimeStamp m_cd_afterpulse_thold{0, 50000};
    TimeStamp m_wp_afterpulse_thold{0, 4000};
    TimeStamp m_cd_last_muon{0, 0};
    TimeStamp m_wp_last_muon{0, 0};

};

#endif // JRAF_MUON_MUONTAGGER_HPP_