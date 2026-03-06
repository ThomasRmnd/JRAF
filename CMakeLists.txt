PKG(AnalysisGroupC

MODULE

FILES

    src/analysis/AccidentalAnalysis.cpp
    src/analysis/Analysis.cpp
    src/analysis/IBDAnalysis.cpp
    src/analysis/MultiplicityAnalysis.cpp
    src/analysis/NeutronAnalysis.cpp

    src/event/Event.cpp
    src/event/EventCache.cpp
    src/event/EventContext.cpp
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