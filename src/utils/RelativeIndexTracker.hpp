#ifndef JRAF_UTILS_RELATIVEINDEXTRACKER_HPP_
#define JRAF_UTILS_RELATIVEINDEXTRACKER_HPP_

#include <filesystem>
#include <string>

#include "SniperKernel/SniperLog.h"

#include "RootIOSvc/RootInputSvc.h"

struct RelativeIndexTracker {

    std::string current;
    std::string next;
    bool change = false;

    std::ptrdiff_t index = -1ll;

    void update(RootInputSvc* iptsvc) {
        ++index;
        std::string filename = std::filesystem::path(iptsvc->getInputStream("EvtNavigator")->streamname()).filename().string();
        if (current.empty()) {
            current = filename;
            next = current;
            index = 0ll;
        }
        if (change) {
            current = next;
            change = false;
            index = 0ll;
        }
        if (filename != current) {
            next = filename;
            change = true;
        }
        LogInfo << "Index of the current file: " << index << '\n';
    }

};

#endif // JRAF_UTILS_RELATIVEINDEXTRACKER_HPP_