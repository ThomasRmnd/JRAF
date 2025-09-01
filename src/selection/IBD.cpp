#include "selection/IBD.hpp"

PromptIBDSelection::PromptIBDSelection(
    const std::shared_ptr<FiducialVolumeSelection>& fiducial,
    const std::shared_ptr<EnergySelection>& energy,
    const std::vector<std::shared_ptr<MuonSelection>>& muons
) :
    m_fiducial{fiducial},
    m_energy{energy},
    m_muons{muons}
{}

bool PromptIBDSelection::isIn(const vertex& vtx) const {
    if (!m_fiducial->isIn(vtx)) return false;
    if (!m_energy->isIn(vtx)) return false;
    for (const std::shared_ptr<MuonSelection>& muon : m_muons) {
        if (muon->isIn(vtx)) return false;
    }
    return true;
}

DelayedIBDSelection::DelayedIBDSelection(
    const std::shared_ptr<FiducialVolumeSelection>& fiducial,
    const std::shared_ptr<EnergySelection>& energy,
    const std::shared_ptr<TimeSelection>& coincidence_time,
    const std::shared_ptr<VolumeSelection>& coincidence_volume,
    const std::vector<std::shared_ptr<MuonSelection>>& muons
) :
    m_fiducial{fiducial},
    m_energy{energy},
    m_coincidence_time{coincidence_time},
    m_coincidence_volume{coincidence_volume},
    m_muons{muons}
{}

bool DelayedIBDSelection::isIn(const vertex& vtx) const {
    if (!m_fiducial->isIn(vtx)) return false;
    if (!m_energy->isIn(vtx)) return false;
    if (!m_coincidence_time->isIn(vtx)) return false;
    if (!m_coincidence_volume->isIn(vtx)) return false;
    for (const std::shared_ptr<MuonSelection>& muon : m_muons) {
        if (muon->isIn(vtx)) return false;
    }
    return true;
}