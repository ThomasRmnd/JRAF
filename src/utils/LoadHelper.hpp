#ifndef RECMUONALG_UTILS_LOADHELPER_HPP_
#define RECMUONALG_UTILS_LOADHELPER_HPP_

#include <list>

#include "Geometry/CdGeom.h"
#include "Geometry/IPMTParamSvc.h"
#include "Event/CalibPmtChannel.h"
#include "RecTools/PmtProp.h"

namespace details {

using CalibPmtList = std::list<JM::CalibPmtChannel*>;

bool loadCol(CalibPmtList::const_iterator clb_first, CalibPmtList::const_iterator clb_last, PmtTable::iterator table_first, PmtTable::iterator table_last);

bool loadCdGeom(PmtTable::iterator first, PmtTable::iterator last, CdGeom* geom, double res_20inch, double res_3inch);

} // namespace details

#endif // RECMUONALG_UTILS_LOADHELPER_HPP_