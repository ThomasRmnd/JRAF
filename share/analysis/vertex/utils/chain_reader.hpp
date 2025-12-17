#ifndef UTILS_CHAIN_HPP_
#define UTILS_CHAIN_HPP_

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>

#include <TChain.h>

class chain_reader {

public:

    bool is_valid() const {
        return m_chain != nullptr && m_chain->GetEntries() > 0;
    }

    std::size_t size() const {
        if (!is_valid()) return 0ul;
        return static_cast<std::size_t>(m_chain->GetEntries());
    }

    bool entry(std::ptrdiff_t n) {
        if (!is_valid()) return false;
        if (n < 0 || static_cast<std::size_t>(n) >= size()) return false;
        m_chain->GetEntry(static_cast<Long64_t>(n));
        return true;
    }

    template<typename _Tp>
    int branch(const std::string& key, _Tp& value) {
        return m_chain->SetBranchAddress(key.c_str(), &value);
    }

private:

    std::shared_ptr<TChain> m_chain;

    chain_reader(const std::string& filepath, const std::string& treename) {
        m_chain = std::make_shared<TChain>(treename.c_str());
        if (m_chain) m_chain->Add(filepath.c_str());
    }

    friend class chain_reader_manager;

};

class chain_reader_manager {

public:

    static std::shared_ptr<chain_reader> retrieve(const std::string& filepath, const std::string& treename) {
        std::string key = filepath + "::" + treename;
        map_type::iterator it = s_cache.find(key);
        if (it != s_cache.end()) {
            if (std::shared_ptr<chain_reader> existing = it->second.lock()) {
                std::cout << "Reusing existing chain for " << key << '\n';
                return existing;
            }
            s_cache.erase(it);
        }
        std::cout << "Creating new chain for " << key << '\n';
        std::shared_ptr<chain_reader> chain(new chain_reader(filepath, treename));
        s_cache[key] = chain;
        return chain;
    }

private:

    typedef std::unordered_map<std::string, std::weak_ptr<chain_reader>> map_type;

    inline static map_type s_cache;

};

#endif // UTILS_CHAIN_HPP_