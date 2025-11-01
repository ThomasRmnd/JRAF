PKG(AnalysisGroupC

MODULE

FILES

    src/analysis/Analysis.cpp
    src/analysis/CdWpCosmoStudy.cpp
    src/analysis/FirstCrossCheckAnalysis.cpp
    src/analysis/IBDWithCylindricalCut.cpp
    src/analysis/IBDWithNeutronVetoStudy.cpp
    src/analysis/MultiplicityWindowCut.cpp
    src/analysis/NavBufferCache.cpp
    src/analysis/NeutronVetoStudy.cpp
    src/analysis/TtCosmoStudy.cpp

    src/event/Event.cpp
    src/event/EventCache.cpp
    src/event/IBD.cpp
    src/event/Track.cpp
    src/event/Vertex.cpp

    src/selection/Energy.cpp
    src/selection/Muon.cpp
    src/selection/Selection.cpp
    src/selection/Time.cpp
    src/selection/Vertex.cpp
    src/selection/Volume.cpp

DEPENDS

    CalibEvent
    RecEvent
    SimEvent
    
    Identifier
    GeometryLib
    RecTools
    UtilsThomas
    BufferMemMgrLib
    Context
    EvtNavigator
    InputReviser
    SNiPER::RootWriterLib

    $<$<BOOL:${WITH_OEC}>:EvtDataTypes>
    $<$<BOOL:${WITH_OEC}>:EvtStore>
    $<$<BOOL:${WITH_OEC}>:OecEvent>

)