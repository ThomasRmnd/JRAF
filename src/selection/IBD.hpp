#ifndef ANALYSISGROUPC_SELECTION_IBD_HPP_
#define ANALYSISGROUPC_SELECTION_IBD_HPP_

#include "selection/Selection.hpp"

#include "selection/Energy.hpp"
#include "selection/Muon.hpp"
#include "selection/Volume.hpp"
#include "selection/Time.hpp"

class IBDSelection : public Selection {

public:

    ~IBDSelection() override = default;

    virtual bool isIn(const vertex& vtx) const = 0;

};

class PromptIBDSelection : public IBDSelection {

public:

    PromptIBDSelection(
        const std::shared_ptr<FiducialVolumeSelection>& fiducial,
        const std::shared_ptr<EnergySelection>& energy,
        const std::vector<std::shared_ptr<MuonSelection>>& muons
    );

    ~PromptIBDSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    std::shared_ptr<FiducialVolumeSelection> m_fiducial;
    std::shared_ptr<EnergySelection> m_energy;
    std::vector<std::shared_ptr<MuonSelection>> m_muons;

};

class DelayedIBDSelection : public IBDSelection {

public:

    DelayedIBDSelection(
        const std::shared_ptr<FiducialVolumeSelection>& fiducial,
        const std::shared_ptr<EnergySelection>& energy,
        const std::shared_ptr<TimeSelection>& coincidence_time,
        const std::shared_ptr<VolumeSelection>& coincidence_volume,
        const std::vector<std::shared_ptr<MuonSelection>>& muons
    );

    ~DelayedIBDSelection() override = default;

    bool isIn(const vertex& vtx) const override;

private:

    std::shared_ptr<FiducialVolumeSelection> m_fiducial;
    std::shared_ptr<EnergySelection> m_energy;
    std::shared_ptr<TimeSelection> m_coincidence_time;
    std::shared_ptr<VolumeSelection> m_coincidence_volume;
    std::vector<std::shared_ptr<MuonSelection>> m_muons;

};

#endif // ANALYSISGROUPC_SELECTION_IBD_HPP_