#ifndef RECMUONALG_LOADER_CD_CDLRANGEFILLER_HPP_
#define RECMUONALG_LOADER_CD_CDLRANGEFILLER_HPP_

#include <list>

#include "Event/CalibPmtChannel.h"

#include "loader/RangeFiller.hpp"

class CdLRangeFiller : public RangeFiller<CdGeom> {

public:

    CdLRangeFiller(const std::string& name, double res_20inch);

    ~CdLRangeFiller() override = default;

    std::size_t size() override;

    bool fill(JM::EvtNavigator* nav) override;

protected:

    double m_res_20inch;

    bool loadGeom() override;
    bool loadCol(const std::list<JM::CalibPmtChannel*>& clb_list);

};

#endif // RECMUONALG_LOADER_CD_CDLRANGEFILLER_HPP_