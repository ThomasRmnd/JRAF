#include "muon/MuonTagger.hpp"

#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(MuonTagger);

MuonTagger::MuonTagger(const std::string& name) :
    ToolBase{name}
{
    declProp("CdMuonTotqThreshold", m_cd_muon_totq_thold = 1000.0);
    declProp("CdOnlyMuonTotqThreshold", m_cd_only_muon_totq_thold = 30000.0);
    declProp("WpMuonTotqThreshold", m_wp_muon_totq_thold = 400.0);
    declProp("WpOnlyMuonTotqThreshold", m_wp_only_muon_totq_thold = 700.0);
    declProp("CdAfterpulseThreshold", m_cd_afterpulse_thold_ns = 50000);
    declProp("WpAfterpulseThreshold", m_wp_afterpulse_thold_ns = 4000);
}

bool MuonTagger::initialize() {
    m_cd_afterpulse_thold = TimeStamp{0, m_cd_afterpulse_thold_ns};
    m_wp_afterpulse_thold = TimeStamp{0, m_wp_afterpulse_thold_ns};
    return true;
}

muon_tagging_result MuonTagger::tag(const muon_tagging_context& ctx) {
    muon_tagging_result result;

    if (
        ctx.calid_cd.totq >= m_cd_muon_totq_thold && 
        ctx.calid_wp.totq >= m_wp_muon_totq_thold && 
        ctx.ts - m_cd_last_muon > m_cd_afterpulse_thold &&
        ctx.ts - m_wp_last_muon > m_wp_afterpulse_thold
    ) {
        m_cd_last_muon = ctx.ts;
        m_wp_last_muon = ctx.ts;
        result.is_possibly_cd_muon = true;
        result.is_possibly_wp_muon = true;
    }
    else if (
        ctx.calid_cd.totq < m_cd_muon_totq_thold && 
        ctx.calid_wp.totq >= m_wp_only_muon_totq_thold &&
        ctx.ts - m_wp_last_muon > m_wp_afterpulse_thold
    ) {
        m_wp_last_muon = ctx.ts;
        result.is_possibly_wp_muon = true;
    }
    else if (
        ctx.calid_cd.totq >= m_cd_only_muon_totq_thold && 
        ctx.calid_wp.totq < m_wp_muon_totq_thold &&
        ctx.ts - (m_cd_last_muon > m_wp_last_muon ? m_cd_last_muon : m_wp_last_muon) > TimeStamp{0, 2000000}
    ) {
        m_cd_last_muon = ctx.ts;
        result.is_possibly_cd_muon = true;
    }

    return result;
}