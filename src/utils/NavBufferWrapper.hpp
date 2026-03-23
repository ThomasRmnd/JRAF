#ifndef JRAF_UTILS_NAVBUFFERWRAPPER_HPP_
#define JRAF_UTILS_NAVBUFFERWRAPPER_HPP_

#include "EvtNavigator/NavBuffer.h"

class NavBufferWrapper : public JM::NavBuffer {

public:

    NavBufferWrapper(JM::NavBuffer& buf) {
        m_dBuf.assign(buf.begin(), buf.end());
        m_iCur = 0;
    }

    ~NavBufferWrapper() override = default;

    void prev() {
        --m_iCur;
    }

    void next() {
        ++m_iCur;
    }

};

#endif // JRAF_UTILS_NAVBUFFERWRAPPER_HPP_