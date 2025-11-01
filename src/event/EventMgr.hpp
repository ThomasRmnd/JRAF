#ifndef ANALYSIS_EVENT_EVENTMGR_HPP_
#define ANALYSIS_EVENT_EVENTMGR_HPP_

#include "EvtNavigator/NavBuffer.h"

class NavBufferWrapper : public JM::NavBuffer {

public:

    NavBufferWrapper(JM::NavBuffer& buf);

    ~NavBufferWrapper() override = default;

    void prev();
    void next();

};

#endif // ANALYSIS_EVENT_EVENTMGR_HPP_