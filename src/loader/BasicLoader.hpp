#ifndef ANALYSISGROUPC_LOADER_BASICLOADER_HPP_
#define ANALYSISGROUPC_LOADER_BASICLOADER_HPP_

#include "loader/Loader.hpp"

class BasicLoader : public Loader {

public:

    using Loader::Loader;

    ~BasicLoader() override = default;

    bool load(JM::NavBuffer* buf) override;

};

#endif // ANALYSISGROUPC_LOADER_BASICLOADER_HPP_