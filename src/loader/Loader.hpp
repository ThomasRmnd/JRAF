#ifndef ANALYSISGROUPC_LOADER_LOADER_HPP_
#define ANALYSISGROUPC_LOADER_LOADER_HPP_

#include <memory>

#include "SniperKernel/ToolBase.h"

#include "EvtNavigator/NavBuffer.h"
#include "Geometry/IRecGeomSvc.hh"
#include "RecTools/PmtProp.h"

#include "loader/RangeFiller.hpp"

/**
 * @class Loader
 * 
 * @brief Base class for PMT table loaders.
 * 
 * It processes `NavBuffer` and fills the corresponding range(s) of PMT using `RangeFiller` for
 * each detector type (CD, WP and TT).
 */
class Loader : public ToolBase {

public:

    Loader(
        const std::string& name, 
        PmtTable* pmt_table, 
        const std::shared_ptr<RangeFiller<CdGeom>>& cd_filler, 
        const std::shared_ptr<RangeFiller<WpGeom>>& wp_filler, 
        const std::shared_ptr<RangeFiller<TtGeom>>& tt_filler, 
        IRecGeomSvc* geom_svc
    );

    virtual ~Loader() = default;

    bool initialize() override;
    bool finalize() override;

    virtual bool load(JM::NavBuffer* buf) = 0;

protected:

    PmtTable* m_pmt_table = nullptr;
    IRecGeomSvc* m_geom_svc = nullptr;

    std::shared_ptr<RangeFiller<CdGeom>> m_cd_filler;
    std::shared_ptr<RangeFiller<WpGeom>> m_wp_filler;
    std::shared_ptr<RangeFiller<TtGeom>> m_tt_filler;

    PmtTable::iterator m_cd_it;
    PmtTable::iterator m_wp_it;
    PmtTable::iterator m_tt_it;

    bool initGeomSvc();
    bool loadNav(JM::EvtNavigator* nav);
    void unload();

};

#endif // ANALYSISGROUPC_LOADER_LOADER_HPP_