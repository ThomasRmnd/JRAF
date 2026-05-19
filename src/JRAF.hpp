#ifndef JRAF_JRAF_HPP_
#define JRAF_JRAF_HPP_

#include "SniperKernel/AlgBase.h"

#include <deque>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <TFile.h>

#include "Context/TimeStamp.h"
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

#include "analysis/Analysis.hpp"
#include "event/Event.hpp"
#include "event/EventBuilder.hpp"
#include "utils/ContextFileTracker.hpp"
#include "utils/DAQTimeSaver.hpp"
#include "utils/FeatureSaver.hpp"
#include "utils/TrackSaver.hpp"
#include "utils/TtRecoFile.hpp"
#include "veto/Veto.hpp"
#include "veto/VetoTimeSaver.hpp"

class JRAF : public AlgBase {

public:

    JRAF(const std::string& name);

    ~JRAF() override = default;

    bool initialize() override;
    bool execute() override;
    bool finalize() override;

private:

    std::size_t m_iEvt = 0ul;
    TimeStamp m_tsEvt{0, 0};
    
    JM::NavBuffer* m_buf;
    RootInputSvc* m_iptSvc;

    std::string m_eventBuilderName;
    EventBuilder* m_eventBuilder;

    ContextFileTracker m_contextTracker;
    BeginningOfJobVetoTracker m_begOfJobVetoTrkr;
    MissingHeaderVetoTracker m_missHdrVetoTrkr;
    BigGapsVetoTracker m_bigGapsVetoTrkr;
    MuonVetoTracker m_muvetoTrkr;

    std::string m_ofilename;
    TFile* m_file;
    DAQTimeSaver m_daqTimeSaver;
    VetoTimeSaver m_vetoTimeSaver;
    std::vector<std::string> m_methods;
    std::vector<std::shared_ptr<Analysis>> m_analyses;

    TrackSaver m_trkSaver;
    FeatureSaver m_featureSaver;

	bool initBufSvc();
    bool initAnalyses();

    // void addFeature(const std::vector<track>& tracks, const TimeStamp& curts, int run_id);

};

#endif // JRAF_JRAF_HPP_