#ifndef RECMUONALG_LOADER_JOINTLOADER_HPP_
#define RECMUONALG_LOADER_JOINTLOADER_HPP_

#include "Context/TimeStamp.h"

#include "loader/Loader.hpp"
#include "utils/DetectorType.hpp"

class JointLoader : public Loader {

public:

    JointLoader(const std::string& name, PmtTable* pmt_table, const std::pair<double, double>& time_window, const std::shared_ptr<RangeFiller<CdGeom>>& cd_loader, const std::shared_ptr<RangeFiller<WpGeom>>& wp_loader, const std::shared_ptr<RangeFiller<TtGeom>>& tt_loader, IRecGeomSvc* geom_svc);

    ~JointLoader() override = default;

    bool load(JM::NavBuffer* buf) override;

protected:

    std::pair<double, double> m_time_window;

    TimeStamp m_cur_ts;
    TimeStamp m_other_ts;

    DetectorType m_cur_type;
    DetectorType m_other_type;

    DetectorType getDetectorType(JM::EvtNavigator* nav);
    void changeRefTime(double diff_ts);
    void changeRefTimeInRange(double diff_ts, PmtTable::iterator first, PmtTable::iterator last);
    void unloadPrev();
    void unloadRange(PmtTable::iterator first, PmtTable::iterator last);

};

#endif // RECMUONALG_LOADER_JOINTLOADER_HPP_