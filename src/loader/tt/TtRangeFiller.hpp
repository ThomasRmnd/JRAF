#ifndef RECMUONALG_LOADER_TT_TTRANGEFILLER_HPP_
#define RECMUONALG_LOADER_TT_TTRANGEFILLER_HPP_

#include <list>

#include "Event/CalibTtChannel.h"
#include "Geometry/ITTGeomSvc.hh"

#include "loader/RangeFiller.hpp"

class TtRangeFiller : public RangeFiller<TtGeom> {

public:

    TtRangeFiller(const std::string& name, double res_tt, ITTGeomSvc* tt_geom);

    ~TtRangeFiller() override = default;

    bool initialize() override;

    std::size_t size() override;

    bool fill(JM::EvtNavigator* nav) override;

protected:

    double m_res_tt;
    ITTGeomSvc* m_tt_geom;

    bool initTTGeomSvc();
    bool loadGeom() override;
    bool loadCol(const std::list<JM::CalibTtChannel*>& clb_list);

};

#endif // RECMUONALG_LOADER_TT_TTRANGEFILLER_HPP_