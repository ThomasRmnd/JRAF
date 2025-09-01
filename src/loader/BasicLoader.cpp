#include "loader/BasicLoader.hpp"

#include "SniperKernel/SniperLog.h"

bool BasicLoader::load(JM::NavBuffer* buf) {
    if (!buf) {
        LogError << "NavBuffer is nullptr\n";
        return false;
    }
    JM::EvtNavigator* nav = buf->curEvt();
    if (!nav) {
        LogError << "Current EvtNavigator is nullptr\n";
        return false;
    }

    unload();
    if (!loadNav(nav)) return false;

    return true;
}