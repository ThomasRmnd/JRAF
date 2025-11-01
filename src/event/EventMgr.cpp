#include "event/EventMgr.hpp"

NavBufferWrapper::NavBufferWrapper(JM::NavBuffer& buf) {
    m_dBuf.assign(buf.begin(), buf.end());
}

void NavBufferWrapper::prev() { 
    --m_iCur; 
}

void NavBufferWrapper::next() { 
    ++m_iCur; 
}