#ifndef RECMUONALG_LOADER_WP_WPRANGEFILLER_HPP_
#define RECMUONALG_LOADER_WP_WPRANGEFILLER_HPP_

#include <list>

#include "Event/CalibPmtChannel.h"

#include "loader/RangeFiller.hpp"

class WpRangeFiller : public RangeFiller<WpGeom> {

public:

    WpRangeFiller(const std::string& name, double res_wp);

    ~WpRangeFiller() override = default;

    std::size_t size() override;

    bool fill(JM::EvtNavigator* nav) override;

protected:

    double m_res_wp;

    bool loadGeom() override;
    bool loadCol(const std::list<JM::CalibPmtChannel*>& clb_list);

};

#endif // RECMUONALG_LOADER_WP_WPRANGEFILLER_HPP_