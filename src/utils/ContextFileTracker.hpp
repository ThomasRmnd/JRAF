#ifndef JRAF_UTILS_CONTEXTFILETRACKER_HPP_
#define JRAF_UTILS_CONTEXTFILETRACKER_HPP_

#include <filesystem>
#include <string>

#include "SniperKernel/SniperLog.h"

#include "RootIOSvc/RootInputSvc.h"

struct ContextFileTracker {

    std::string prevctx;
    std::string nextctx;

    std::string current;
    std::string next;
    bool change = false;

    bool isTarget(RootInputSvc* iptSvc) {
        if (current.empty()) {
            current = std::filesystem::path(iptSvc->getInputStream("EvtNavigator")->streamname()).filename().string();
            next = current;
        }
        if (change) {
            current = next;
            change = false;
        }
        std::string filename = std::filesystem::path(iptSvc->getInputStream("EvtNavigator")->streamname()).filename().string();
        if (filename != current) {
            next = filename;
            change = true;
        }
        LogInfo << "Current filename: " << current << '\n';
        LogInfo << "Previous filename: " << prevctx << '\n';
        LogInfo << "Next filename: " << nextctx << '\n';
        return (current != prevctx) && (current != nextctx);
    }

};

#endif // JRAF_UTILS_CONTEXTFILETRACKER_HPP_