PKG(AnalysisGroupC

MODULE

FILES

    src/analysis/Analysis.cpp
    src/analysis/CdWpCosmoStudy.cpp
    src/analysis/FirstCrossCheckAnalysis.cpp
    src/analysis/IBDWithCylindricalCut.cpp
    src/analysis/IBDWithNeutronVetoStudy.cpp
    src/analysis/MultiplicityWindowCut.cpp
    src/analysis/NeutronVetoStudy.cpp
    src/analysis/TtCosmoStudy.cpp

    src/event/Event.cpp
    src/event/EventCache.cpp
    src/event/IBD.cpp
    src/event/Track.cpp
    src/event/Vertex.cpp

    src/loader/BasicLoader.cpp
    src/loader/Loader.cpp
    src/loader/JointLoader.cpp
    src/loader/cd/CdLRangeFiller.cpp
    src/loader/cd/CdRangeFiller.cpp
    src/loader/cd/CdSRangeFiller.cpp
    src/loader/tt/TtRangeFiller.cpp
    src/loader/wp/AkiraWpRangeFiller.cpp
    src/loader/wp/WpRangeFiller.cpp

    src/selection/Energy.cpp
    src/selection/IBD.cpp
    src/selection/Muon.cpp
    src/selection/Selection.cpp
    src/selection/Time.cpp
    src/selection/Volume.cpp

    src/utils/LoadHelper.cpp

DEPENDS

    CalibEvent
    RecEvent
    SimEvent
    
    Identifier
    GeometryLib
    RecTools
    BufferMemMgrLib
    Context
    EvtNavigator
    InputReviser
    SNiPER::RootWriterLib

    $<$<BOOL:${WITH_OEC}>:EvtDataTypes>
    $<$<BOOL:${WITH_OEC}>:EvtStore>
    $<$<BOOL:${WITH_OEC}>:OecEvent>

)