#include "utils/LoadHelper.hpp"

#include "SniperKernel/SniperLog.h"

#include "Identifier/CdID.h"
#include "Identifier/Identifier.h"
#include "Identifier/JunoDetectorID.h"
#include "Identifier/TtID.h"
#include "Identifier/WpID.h"

bool details::loadCol(CalibPmtList::const_iterator clb_first, CalibPmtList::const_iterator clb_last, PmtTable::iterator table_first, PmtTable::iterator table_last) {
    for (CalibPmtList::const_iterator it = clb_first; it != clb_last; ++it) {
        Identifier id = Identifier((*it)->pmtId());
        int i;
        
        if (JunoDetectorID::isCd(id)) {
            i = CdID::module(id);
        }
        else if (JunoDetectorID::isWp(id)) {
            i = WpID::module(id);
        }
        else {
            LogError << "PMT (id: " << id << ") is not a CD nor a WP PMT\n";
            return false;
        }

        if (i < 0 || std::distance(table_first, table_last) <= i) {
            LogWarn << "PMT (id: " << id << ") is out of range. Index is " << i << ", while the range is " << std::distance(table_first, table_last) << '\n';
            continue;
        }

        PmtTable::iterator it_table = table_first + static_cast<std::ptrdiff_t>(i);
        it_table->used = true;
        it_table->q = (*it)->nPE();
        it_table->fht = (*it)->firstHitTime();
        it_table->hitq = (*it)->charge();
        it_table->hittime = (*it)->time();
    }
    return true;
}

bool details::loadCdGeom(PmtTable::iterator first, PmtTable::iterator last, CdGeom* geom, double res_20inch, double res_3inch) {
    for (PmtTable::iterator it = first; it != last; ++it) {
        Identifier id = Identifier(CdID::id(CdID::offset_cd_spmt(std::distance(first, it)), 0));
        PmtGeom* pmt = geom->getPmt(id);
        if (!pmt) {
            LogError << "Failed to get PMT (id: " << id << ") for CD\n";
            return false;
        }
        it->pmtid = id.getValue();
        it->value = id.getValue();
        it->used = false;
        it->pos = pmt->getCenter();
        it->loc = 1;
        if (CdID::is20inch(id)) {
            it->res = res_20inch;
            it->type = PmtType::PMT_20INCH;
        } 
        else if (CdID::is3inch(id)) {
            it->res = res_3inch;
            it->type = PmtType::PMT_3INCH;
        }
        else {
            LogError << "PMT (id: " << id << ") is not a 20-inch nor a 3-inch PMT\n";
            return false;
        }
    }
    return true;
}