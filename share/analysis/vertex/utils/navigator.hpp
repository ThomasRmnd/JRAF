#ifndef UTILS_NAVIGATOR_HPP_
#define UTILS_NAVIGATOR_HPP_

#include <iostream>
#include <memory>
#include <string>
#include <typeinfo>
#include <unordered_map>

#include "utils/chain_reader.hpp"
#include "utils/event.hpp"
#include "utils/vector_reader.hpp"

class navigator_base {

public:

    navigator_base(const std::string& filepath, const std::string& treename) {
        m_chain = chain_reader_manager::retrieve(filepath, treename);
        if (!m_chain->is_valid()) {
            std::cerr << "Chain of filepath " << filepath << " and treename " << treename << " is not valid\n";
            return;
        }

        m_chain->branch("run_id", run_id);

        m_chain->branch("posx_p", prompt.pos.x);
        m_chain->branch("posy_p", prompt.pos.y);
        m_chain->branch("posz_p", prompt.pos.z);
        m_chain->branch("e_p", prompt.e);
        m_chain->branch("sec_p", prompt.ts.sec);
        m_chain->branch("nsec_p", prompt.ts.nsec);

        m_chain->branch("totq_p", meta_prompt.totq);
        m_chain->branch("meanq_p", meta_prompt.meanq);
        m_chain->branch("stdq_p", meta_prompt.stdq);
        m_chain->branch("minq_p", meta_prompt.minq);
        m_chain->branch("maxq_p", meta_prompt.maxq);
        m_chain->branch("meant_p", meta_prompt.meant);
        m_chain->branch("stdt_p", meta_prompt.stdt);
        m_chain->branch("npmt_p", meta_prompt.npmt);
        m_chain->branch("nhit_p", meta_prompt.nhit);
        m_chain->branch("meanhit_p", meta_prompt.meanhit);
        m_chain->branch("stdhit_p", meta_prompt.stdhit);

        m_chain->branch("posx_d", delayed.pos.x);
        m_chain->branch("posy_d", delayed.pos.y);
        m_chain->branch("posz_d", delayed.pos.z);
        m_chain->branch("e_d", delayed.e);
        m_chain->branch("sec_d", delayed.ts.sec);
        m_chain->branch("nsec_d", delayed.ts.nsec);

        m_chain->branch("totq_d", meta_delayed.totq);
        m_chain->branch("meanq_d", meta_delayed.meanq);
        m_chain->branch("stdq_d", meta_delayed.stdq);
        m_chain->branch("minq_d", meta_delayed.minq);
        m_chain->branch("maxq_d", meta_delayed.maxq);
        m_chain->branch("meant_d", meta_delayed.meant);
        m_chain->branch("stdt_d", meta_delayed.stdt);
        m_chain->branch("npmt_d", meta_delayed.npmt);
        m_chain->branch("nhit_d", meta_delayed.nhit);
        m_chain->branch("meanhit_d", meta_delayed.meanhit);
        m_chain->branch("stdhit_d", meta_delayed.stdhit);
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

        m_chain->branch("posx_n", posx_n.data());
        m_chain->branch("posy_n", posy_n.data());
        m_chain->branch("posz_n", posz_n.data());
        m_chain->branch("e_n", e_n.data());
        m_chain->branch("sec_n", sec_n.data());
        m_chain->branch("nsec_n", nsec_n.data());

        m_chain->branch("totq_n", totq_n.data());
        m_chain->branch("meanq_n", meanq_n.data());
        m_chain->branch("stdq_n", stdq_n.data());
        m_chain->branch("minq_n", minq_n.data());
        m_chain->branch("maxq_n", maxq_n.data());
        m_chain->branch("meant_n", meant_n.data());
        m_chain->branch("stdt_n", stdt_n.data());
        m_chain->branch("npmt_n", npmt_n.data());
        m_chain->branch("nhit_n", nhit_n.data());
        m_chain->branch("meanhit_n", meanhit_n.data());
        m_chain->branch("stdhit_n", stdhit_n.data());

        m_chain->branch("posx_mult", posx_mult.data());
        m_chain->branch("posy_mult", posy_mult.data());
        m_chain->branch("posz_mult", posz_mult.data());
        m_chain->branch("e_mult", e_mult.data());
        m_chain->branch("sec_mult", sec_mult.data());
        m_chain->branch("nsec_mult", nsec_mult.data());
        m_chain->branch("mult_type", mult_type.data());

        m_chain->branch("totq_mult", totq_mult.data());
        m_chain->branch("meanq_mult", meanq_mult.data());
        m_chain->branch("stdq_mult", stdq_mult.data());
        m_chain->branch("minq_mult", minq_mult.data());
        m_chain->branch("maxq_mult", maxq_mult.data());
        m_chain->branch("meant_mult", meant_mult.data());
        m_chain->branch("stdt_mult", stdt_mult.data());
        m_chain->branch("npmt_mult", npmt_mult.data());
        m_chain->branch("nhit_mult", nhit_mult.data());
        m_chain->branch("meanhit_mult", meanhit_mult.data());
        m_chain->branch("stdhit_mult", stdhit_mult.data());

        m_chain->branch("method_mu", method_mu.data());
        m_chain->branch("loc_mu", loc_mu.data());
        m_chain->branch("posx_mu", posx_mu.data());
        m_chain->branch("posy_mu", posy_mu.data());
        m_chain->branch("posz_mu", posz_mu.data());
        m_chain->branch("dirx_mu", dirx_mu.data());
        m_chain->branch("diry_mu", diry_mu.data());
        m_chain->branch("dirz_mu", dirz_mu.data());
        m_chain->branch("totq_cd_mu", totq_cd_mu.data());
        m_chain->branch("totq_wp_mu", totq_wp_mu.data());
        m_chain->branch("sec_mu", sec_mu.data());
        m_chain->branch("nsec_mu", nsec_mu.data());
        m_chain->branch("quality_mu", quality_mu.data());
    }

    virtual ~basic_navigator() override = default;

    vector_reader<double> posx_n;
    vector_reader<double> posy_n;
    vector_reader<double> posz_n;
    vector_reader<double> e_n;
    vector_reader<time_t> sec_n;
    vector_reader<int> nsec_n;

    vector_reader<double> totq_n;
    vector_reader<double> meanq_n;
    vector_reader<double> stdq_n;
    vector_reader<double> minq_n;
    vector_reader<double> maxq_n;
    vector_reader<double> meant_n;
    vector_reader<double> stdt_n;
    vector_reader<std::size_t> npmt_n;
    vector_reader<std::size_t> nhit_n;
    vector_reader<double> meanhit_n;
    vector_reader<double> stdhit_n;

    vector_reader<double> posx_mult;
    vector_reader<double> posy_mult;
    vector_reader<double> posz_mult;
    vector_reader<double> e_mult;
    vector_reader<time_t> sec_mult;
    vector_reader<int> nsec_mult;
    vector_reader<int> mult_type;

    vector_reader<double> totq_mult;
    vector_reader<double> meanq_mult;
    vector_reader<double> stdq_mult;
    vector_reader<double> minq_mult;
    vector_reader<double> maxq_mult;
    vector_reader<double> meant_mult;
    vector_reader<double> stdt_mult;
    vector_reader<std::size_t> npmt_mult;
    vector_reader<std::size_t> nhit_mult;
    vector_reader<double> meanhit_mult;
    vector_reader<double> stdhit_mult;

    vector_reader<std::string> method_mu;
    vector_reader<int> loc_mu;
    vector_reader<double> posx_mu;
    vector_reader<double> posy_mu;
    vector_reader<double> posz_mu;
    vector_reader<double> dirx_mu;
    vector_reader<double> diry_mu;
    vector_reader<double> dirz_mu;
    vector_reader<double> totq_cd_mu;
    vector_reader<double> totq_wp_mu;
    vector_reader<time_t> sec_mu;
    vector_reader<int> nsec_mu;
    vector_reader<double> quality_mu;

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