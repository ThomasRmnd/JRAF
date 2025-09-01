#ifndef ANALYSISGROUPC_LOADER_RANGEFILLER_HPP_
#define ANALYSISGROUPC_LOADER_RANGEFILLER_HPP_

#include "SniperKernel/ToolBase.h"

#include "EvtNavigator/NavBuffer.h"
#include "Geometry/CdGeom.h"
#include "Geometry/TtGeom.h"
#include "Geometry/WpGeom.h"
#include "RecTools/PmtProp.h"

/**
 * @class RangeFiller
 * @tparam GeomType The geometry type, which defines the detector type.
 * 
 * @brief Base class used to manage a range of the PMT table. 
 * 
 * To correctly initialize this class:
 * 1. Set the geometry with `setGeom`.
 * 2. Get the number of PMTs with `size`.
 * 3. Set the range with `setRange`.
 * 4. Initialize with `initialize`.
 * 5. Use fill for each event
 * 
 * It processes `EvtNavigator` to fill the range of PMT.
 */
template <typename GeomType>
class RangeFiller : public ToolBase {

public:

    using ToolBase::ToolBase;

    virtual ~RangeFiller() = default;

    virtual bool initialize() override {
        if (!initGeomSvc()) return false;
        if (!checkRange()) return false;
        if (!loadGeom()) return false;
        return true;
    }
    
    virtual bool finalize() override { return true; }

    // Set the geometry.
    void setGeom(GeomType* geom) { m_geom = geom; }
    
    // Get the size of the PMT table. This method need to be called after setGeom.
    virtual std::size_t size() = 0;
    
    // Set the range.
    void setRange(PmtTable::iterator first, PmtTable::iterator last) {
        m_first = first;
        m_last = last;
        m_range_set = true;
    }

    virtual bool fill(JM::EvtNavigator* nav) = 0;

protected:

    PmtTable::iterator m_first;
    PmtTable::iterator m_last;
    bool m_range_set = false;
    GeomType* m_geom;

    bool initGeomSvc() {
        if (!m_geom) {
            LogError << "Geom has not been set." << std::endl;
            return false;
        }
        return true;
    }

    bool checkRange() {
        if (!m_range_set) {
            LogError << "Range has not been set." << std::endl;
            return false;
        }
        return true;
    }

    virtual bool loadGeom() = 0;

};


#endif // ANALYSISGROUPC_LOADER_RANGEFILLER_HPP_