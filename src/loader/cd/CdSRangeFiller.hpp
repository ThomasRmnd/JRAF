#ifndef RECMUONALG_LOADER_CD_CDSRANGEFILLER_HPP_
#define RECMUONALG_LOADER_CD_CDSRANGEFILLER_HPP_

#include <list>

#include "Event/CalibPmtChannel.h"

#include "loader/RangeFiller.hpp"

class CdSRangeFiller : public RangeFiller<CdGeom> {

public:

    CdSRangeFiller(const std::string& name, double res_3inch);

    ~CdSRangeFiller() override = default;

    std::size_t size() override;

    bool fill(JM::EvtNavigator* nav) override;

protected:

    double m_res_3inch;
    
    bool loadGeom() override;
    bool loadCol(const std::list<JM::CalibPmtChannel*>& clb_list);

};

#endif // RECMUONALG_LOADER_CD_CDSRANGEFILLER_HPP_