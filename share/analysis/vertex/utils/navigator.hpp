#ifndef UTILS_NAVIGATOR_HPP_
#define UTILS_NAVIGATOR_HPP_

#include <iostream>
#include <memory>
#include <string>
#include <typeinfo>
#include <unordered_map>

#include "utils/chain_reader.hpp"
#include "utils/event.hpp"

class navigator_base {

public:

    navigator_base(const std::string& filepath, const std::string& treename) {
        m_chain = chain_reader_manager::retrieve(filepath, treename);
        if (!m_chain->is_valid()) {
            std::cerr << "Chain of filepath " << filepath << " and treename " << treename << " is not valid\n";
            return;
        }

        // m_chain->branch("run_id", run_id);

        m_chain->branch("posx_p", prompt.pos.x);
        m_chain->branch("posy_p", prompt.pos.y);
        m_chain->branch("posz_p", prompt.pos.z);
        m_chain->branch("e_p", prompt.e);
        m_chain->branch("totq_p", prompt.q);
        m_chain->branch("sec_p", prompt.ts.sec);
        m_chain->branch("nsec_p", prompt.ts.nsec);

        m_chain->branch("meanq_p", meta_prompt.meanq);
        m_chain->branch("stdq_p", meta_prompt.stdq);
        m_chain->branch("minq_p", meta_prompt.minq);
        m_chain->branch("maxq_p", meta_prompt.maxq);
        m_chain->branch("nhit_p", meta_prompt.nhit);
        m_chain->branch("meant_p", meta_prompt.meant);
        m_chain->branch("stdt_p", meta_prompt.stdt);

        m_chain->branch("posx_d", delayed.pos.x);
        m_chain->branch("posy_d", delayed.pos.y);
        m_chain->branch("posz_d", delayed.pos.z);
        m_chain->branch("e_d", delayed.e);
        m_chain->branch("totq_d", delayed.q);
        m_chain->branch("sec_d", delayed.ts.sec);
        m_chain->branch("nsec_d", delayed.ts.nsec);

        m_chain->branch("meanq_d", meta_delayed.meanq);
        m_chain->branch("stdq_d", meta_delayed.stdq);
        m_chain->branch("minq_d", meta_delayed.minq);
        m_chain->branch("maxq_d", meta_delayed.maxq);
        m_chain->branch("nhit_d", meta_delayed.nhit);
        m_chain->branch("meant_d", meta_delayed.meant);
        m_chain->branch("stdt_d", meta_delayed.stdt);
    }

    virtual ~navigator_base() = default;

    bool is_valid() const { return m_chain->is_valid(); }

    std::size_t size() const { return m_chain->size(); }

    virtual bool entry(std::ptrdiff_t n) { return m_chain->entry(n); }

    int run_id;

    vertex prompt;
    vertex delayed;

    vertex_metadata meta_prompt;
    vertex_metadata meta_delayed;

protected:

    std::shared_ptr<chain_reader> m_chain;

};

class basic_navigator : public navigator_base {

public:

    basic_navigator(const std::string& filepath, const std::string& treename) : 
        navigator_base(filepath, treename) 
    {
        if (!is_valid()) return;

        m_chain->branch("posx_n", posx_n);
        m_chain->branch("posy_n", posy_n);
        m_chain->branch("posz_n", posz_n);
        m_chain->branch("e_n", e_n);
        m_chain->branch("totq_n", totq_n);
        m_chain->branch("sec_n", sec_n);
        m_chain->branch("nsec_n", nsec_n);

        m_chain->branch("posx_mult", posx_mult);
        m_chain->branch("posy_mult", posy_mult);
        m_chain->branch("posz_mult", posz_mult);
        m_chain->branch("e_mult", e_mult);
        m_chain->branch("totq_mult", totq_mult);
        m_chain->branch("sec_mult", sec_mult);
        m_chain->branch("nsec_mult", nsec_mult);
        m_chain->branch("mult_type", mult_type);

        m_chain->branch("method_mu", method_mu);
        m_chain->branch("loc_mu", loc_mu);
        m_chain->branch("posx_mu", posx_mu);
        m_chain->branch("posy_mu", posy_mu);
        m_chain->branch("posz_mu", posz_mu);
        m_chain->branch("dirx_mu", dirx_mu);
        m_chain->branch("diry_mu", diry_mu);
        m_chain->branch("dirz_mu", dirz_mu);
        m_chain->branch("totq_mu", totq_mu);
        m_chain->branch("sec_mu", sec_mu);
        m_chain->branch("nsec_mu", nsec_mu);
        m_chain->branch("quality_mu", quality_mu);
    }

    virtual ~basic_navigator() override = default;

    std::vector<double>* posx_n = nullptr;
    std::vector<double>* posy_n = nullptr;
    std::vector<double>* posz_n = nullptr;
    std::vector<double>* e_n = nullptr;
    std::vector<double>* totq_n = nullptr;
    std::vector<time_t>* sec_n = nullptr;
    std::vector<int>* nsec_n = nullptr;

    std::vector<double>* posx_mult = nullptr;
    std::vector<double>* posy_mult = nullptr;
    std::vector<double>* posz_mult = nullptr;
    std::vector<double>* e_mult = nullptr;
    std::vector<double>* totq_mult = nullptr;
    std::vector<time_t>* sec_mult = nullptr;
    std::vector<int>* nsec_mult = nullptr;
    std::vector<int>* mult_type = nullptr;

    std::vector<std::string>* method_mu = nullptr;
    std::vector<int>* loc_mu = nullptr;
    std::vector<double>* posx_mu = nullptr;
    std::vector<double>* posy_mu = nullptr;
    std::vector<double>* posz_mu = nullptr;
    std::vector<double>* dirx_mu = nullptr;
    std::vector<double>* diry_mu = nullptr;
    std::vector<double>* dirz_mu = nullptr;
    std::vector<double>* totq_mu = nullptr;
    std::vector<time_t>* sec_mu = nullptr;
    std::vector<int>* nsec_mu = nullptr;
    std::vector<double>* quality_mu = nullptr;

};

using ibd_navigator = basic_navigator;

class navigator_manager {

public:

    template<typename _Nav>
    static std::shared_ptr<_Nav> retrieve(const std::string& filepath, const std::string& treename) {
        std::string key = typeid(_Nav).name();
        key += "::" + filepath + "::" + treename;

        map_type::iterator it = s_cache.find(key);
        if (it != s_cache.end()) {
            if (std::shared_ptr<navigator_base> existing = it->second.lock()) {
                if (std::shared_ptr<_Nav> typed_ptr = std::dynamic_pointer_cast<_Nav>(existing)) {
                    std::cout << "Reusing navigator " << key << '\n';
                    return typed_ptr;
                }
            }
            s_cache.erase(it);
        }

        std::cout << "Creating new Navigator for " << key << '\n';
        std::shared_ptr<_Nav> nav(new _Nav(filepath, treename));
        
        if (!nav->is_valid()) {
            std::cerr << "Failed to create and validate Navigator " << key << '\n';
            return nullptr;
        }

        s_cache[key] = nav; 
        return nav;
    }

private:

    typedef std::unordered_map<std::string, std::weak_ptr<navigator_base>> map_type;
    inline static map_type s_cache;

    navigator_manager() = delete;

};

#endif // UTILS_NAVIGATOR_HPP_