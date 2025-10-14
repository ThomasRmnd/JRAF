#ifndef RECMUONALG_LOADER_WP_AKIRAWPRANGEFILLER_HPP_
#define RECMUONALG_LOADER_WP_AKIRAWPRANGEFILLER_HPP_

#include "loader/RangeFiller.hpp"

#include <list>

#include "Event/CalibPmtChannel.h"

class AkiraWpRangeFiller : public RangeFiller<WpGeom> {

public:

    AkiraWpRangeFiller(const std::string& name);
    AkiraWpRangeFiller(const std::string& name, double res_wp);

    ~AkiraWpRangeFiller() override = default;

    std::size_t size() override;

    bool fill(JM::EvtNavigator* nav) override;

protected:

    double m_res_wp;

    bool loadGeom() override;
    bool loadCol(const std::list<JM::CalibPmtChannel*>& clb_list);

};

#endif // RECMUONALG_LOADER_WP_AKIRAWPRANGEFILLER_HPP_