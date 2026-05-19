#ifndef JRAF_EVENT_EVENTBUILDER_HPP_
#define JRAF_EVENT_EVENTBUILDER_HPP_

#include <list>

#include "SniperKernel/AlgBase.h"
#include "SniperKernel/ToolBase.h"

#include "Event/CalibPmtChannel.h"
#include "Event/CdTrackRecHeader.h"
#include "Event/CdVertexRecHeader.h"
#include "Event/OecHeader.h"
#include "Event/TtRecHeader.h"
#include "Event/WpRecHeader.h"
#include "EvtNavigator/NavBuffer.h"
#include "RecTools/IRecMuonTool.h"
#include "UtilsThomas/loader/Loader.hpp"
#include "UtilsThomas/utils/DetectorType.hpp"

#include "event/CalibrationContext.hpp"
#include "event/Track.hpp"
#include "event/Vertex.hpp"
#include "muon/MuonTagger.hpp"
#include "utils/TtRecoFile.hpp"

class EventBuilder : public ToolBase {

public:

    EventBuilder(const std::string& name);

    virtual ~EventBuilder() override = default;

    bool initialize() override;
    bool finalize() override;

    bool configure(AlgBase* alg);
    virtual bool build(JM::NavBuffer* buf);

protected:

    std::string m_mutagger_name;
    MuonTagger* m_mutagger;

    std::string m_loader_name;
    std::string m_cdfiller_name;
    std::string m_wpfiller_name;
    std::string m_ttfiller_name;
    Loader* m_loader;
    PmtTable m_pmt_table;
    Params m_params;

    std::string m_rectool_name;
    IRecMuonTool* m_rectool;
    TtRecoFile m_tt_reco_file;

    DetectorType getDetectorType(JM::EvtNavigator* nav);
    calibration_context getCalibrationContext(const std::list<JM::CalibPmtChannel*>& clb_list);
    
    void addTrack(RecTrks& rec_tracks, const std::string& method, double totq_cd, double totq_wp, const TimeStamp& ts, const track::loc& det, std::vector<track>& tracks);
    void addTrack(JM::CdTrackRecHeader* cdt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks);
    void addTrack(JM::WpRecHeader* wpt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks);
    void addTrack(JM::TtRecHeader* ttt_hdr, const std::string& method, const TimeStamp& ts, std::vector<track>& tracks);
    
    void addVertex(JM::OecHeader* oec_hdr, const std::string& method, const TimeStamp& ts, const calibration_context& calib, std::vector<vertex>& vertices);
    void addVertex(JM::CdVertexRecHeader* cdv_hdr, const std::string& method, const TimeStamp& ts, const calibration_context& calib, std::vector<vertex>& vertices);

    int getTtLayerId(double z);
    void addTtToTrack(std::vector<track>& tracks, const TimeStamp& curts);

};

#endif // JRAF_EVENT_EVENTBUILDER_HPP_