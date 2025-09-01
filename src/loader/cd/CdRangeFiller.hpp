#ifndef ANALYSISGROUPC_LOADER_CD_CDRANGEFILLER_HPP_
#define ANALYSISGROUPC_LOADER_CD_CDRANGEFILLER_HPP_

#include <list>

#include "Event/CalibPmtChannel.h"

#include "loader/RangeFiller.hpp"

class CdRangeFiller : public RangeFiller<CdGeom> {

public:

    CdRangeFiller(const std::string& name, double res_20inch, double res_3inch);

    ~CdRangeFiller() override = default;

    std::size_t size() override;

    bool fill(JM::EvtNavigator* nav) override;

protected:

    double m_res_20inch;
    double m_res_3inch;

    bool loadGeom() override;
    bool loadCol(const std::list<JM::CalibPmtChannel*>& clb_list);

};

#endif // ANALYSISGROUPC_LOADER_CD_CDRANGEFILLER_HPP_