#include "loader/JointLoader.hpp"

#include <bitset>
#include <numeric>

#include "Event/CdLpmtCalibHeader.h"
#include "Event/CdSpmtCalibHeader.h"
#include "Event/TtCalibHeader.h"
#include "Event/WpCalibHeader.h"
#include "EvtNavigator/EvtNavHelper.h"
#include "SniperKernel/ToolFactory.h"

JointLoader::JointLoader(const std::string& name, PmtTable* pmt_table, const std::pair<double, double>& time_window, const std::shared_ptr<RangeFiller<CdGeom>>& cd_loader, const std::shared_ptr<RangeFiller<WpGeom>>& wp_loader, const std::shared_ptr<RangeFiller<TtGeom>>& tt_loader, IRecGeomSvc* geom_svc) :
    Loader{name, pmt_table, cd_loader, wp_loader, tt_loader, geom_svc},
    m_time_window{time_window},
    m_cur_ts{},
    m_other_ts{},
    m_cur_type{DetectorType::UNKNOWN},
    m_other_type{DetectorType::UNKNOWN}
{}

bool JointLoader::load(JM::NavBuffer* buf) {
    if (!buf) {
        LogError << "NavBuffer is nullptr\n";
        return false;
    }
    JM::EvtNavigator* nav = buf->curEvt();
    if (!nav) {
        LogError << "Current EvtNavigator is nullptr\n";
        return false;
    }

    unloadPrev();

    m_cur_ts = TimeStamp{nav->TimeStamp().GetTimeSpec()};
    m_cur_type = getDetectorType(nav);
    if (m_cur_type == DetectorType::UNKNOWN) {
        LogWarn << "Unknown entry type at " << m_cur_ts << '\n';
        return true;
    }

    if (!loadNav(nav)) return false;

    for (JM::NavBuffer::Iterator it = buf->begin(); it != buf->end(); ++it) {
        if (it == buf->current()) continue;
        nav = it->get();
        if (!nav) {
            LogError << "EvtNavigator is nullptr\n";
            return false;
        }
        m_other_ts = TimeStamp{nav->TimeStamp().GetTimeSpec()};
        m_other_type = getDetectorType(nav);
        if (m_other_type == DetectorType::UNKNOWN) {
            LogWarn << "Unkown entry type at " << m_other_ts << '\n';
            continue;
        }
        LogDebug << "Current entry: (type: " << static_cast<int>(m_cur_type) << ", ts: " << m_cur_ts << ")\n";
        LogDebug << "Other entry: (type: " << static_cast<int>(m_other_type) << ", ts: " << m_other_ts << ")\n";
        if ( (m_cur_type & m_other_type) != DetectorType::UNKNOWN ) continue; // current event has already assotiated an event with same type
        TimeStamp diff_ts = m_cur_ts - m_other_ts;
        double diff_ts_ns = static_cast<double>(diff_ts.GetSec()) * 1.0e9 + static_cast<double>(diff_ts.GetNanoSec());
        if (diff_ts_ns < m_time_window.first || m_time_window.second < diff_ts_ns) continue;
        LogInfo << "Entries are associated: " << m_cur_ts << " - " << m_other_ts << " = " << diff_ts_ns << " ns\n";

        if (!loadNav(nav)) return false;

        changeRefTime(diff_ts_ns);
        m_cur_type |= m_other_type;
    }

    return true;
}

DetectorType JointLoader::getDetectorType(JM::EvtNavigator* nav) {
    DetectorType type = DetectorType::UNKNOWN;

    JM::CdLpmtCalibHeader* cdl_hdr = JM::getHeaderObject<JM::CdLpmtCalibHeader>(nav);
    if (cdl_hdr && cdl_hdr->event()) type |= DetectorType::CD;
    JM::CdSpmtCalibHeader* cds_hdr = JM::getHeaderObject<JM::CdSpmtCalibHeader>(nav);
    if (cds_hdr && cds_hdr->event()) type |= DetectorType::CD;
    JM::WpCalibHeader* wp_hdr = JM::getHeaderObject<JM::WpCalibHeader>(nav);
    if (wp_hdr && wp_hdr->event()) type |= DetectorType::WP;
    JM::TtCalibHeader* tt_hdr = JM::getHeaderObject<JM::TtCalibHeader>(nav);
    if (tt_hdr && tt_hdr->event()) type |= DetectorType::TT;

    return type;
}

void JointLoader::changeRefTime(double diff_ts) {
    if ( (m_other_type & DetectorType::CD) == DetectorType::CD ) changeRefTimeInRange(diff_ts, m_cd_it, m_wp_it);
    if ( (m_other_type & DetectorType::WP) == DetectorType::WP ) changeRefTimeInRange(diff_ts, m_wp_it, m_tt_it);
    if ( (m_other_type & DetectorType::TT) == DetectorType::TT ) changeRefTimeInRange(diff_ts, m_tt_it, m_pmt_table->end());
}

void JointLoader::changeRefTimeInRange(double diff_ts, PmtTable::iterator first, PmtTable::iterator last) {
    float diff_ts_f32 = static_cast<float>(diff_ts);
    for (PmtTable::iterator it = first; it != last; ++it) {
        if (!it->used) continue;
        it->fht -= diff_ts;
        for (float& time : it->hittime) {
            time -= diff_ts_f32;
        }
    }
}

void JointLoader::unloadPrev() {
    if ( (m_cur_type & DetectorType::CD) == DetectorType::CD ) unloadRange(m_cd_it, m_wp_it);
    if ( (m_cur_type & DetectorType::WP) == DetectorType::WP ) unloadRange(m_wp_it, m_tt_it);
    if ( (m_cur_type & DetectorType::TT) == DetectorType::TT ) unloadRange(m_tt_it, m_pmt_table->end());
}

void JointLoader::unloadRange(PmtTable::iterator first, PmtTable::iterator last) {
    for (PmtTable::iterator it = first; it != last; ++it) {
        it->used = false;
    }
}