#ifndef ANALYSIS_ANALYSIS_HPP_
#define ANALYSIS_ANALYSIS_HPP_

#include "utils/navigator.hpp"

class analysis_base {

public:

    analysis_base(const std::string& name) :
        m_name(name)
    {}

    virtual ~analysis_base() = default;

    virtual std::shared_ptr<navigator_base> navigator() const = 0;

    virtual bool selection() = 0;
    virtual bool process() = 0;
    
    virtual bool save() = 0;
    virtual void result() = 0;

protected:

    std::string m_name;

};

class analysis_registry {

public:

    typedef std::vector<std::shared_ptr<analysis_base>> vector_type;
    typedef std::map<std::shared_ptr<navigator_base>, vector_type> map_type;

    bool book(const std::shared_ptr<analysis_base>& analysis) {
        if (!analysis) {
            std::cerr << "Cannot register analysis\n";
            return false;
        }
        m_registry[analysis->navigator()].push_back(analysis);
        return true;
    }

private:

    map_type m_registry;

    friend class analysis_manager;

};

class analysis_manager {

public:

    analysis_manager(analysis_registry& reg) :
        m_reg{reg}
    {}

    bool run() {
        if (m_reg.m_registry.empty()) {
            std::cout << "WARNING: Analysis registry is empty. Exiting run\n";
            return true;
        }
        
        std::cout << "\n--- Starting Analysis Loop Over " << m_reg.m_registry.size() << " Data Groups ---\n";
        for (auto const& [nav, analyses] : m_reg.m_registry) {
            if (!nav->is_valid()) {
                std::cerr << "\nWARNING: Navigator for this group is invalid. Skipping group\n";
                continue;
            }

            Long64_t entries = nav->size();
            std::cout << "\n[Group Start] Processing " << analyses.size() << " analyses over " << entries << " entries\n";
            
            for (Long64_t k = 0; k < entries; ++k) {
                nav->entry(k); 

                for (const auto& analysis : analyses) {
                    if (analysis->selection()) {
                        if (!analysis->process()) return false;
                    }
                }

                if ((k + 1) % 1000 == 0) {
                     std::cout << "  [Group Status] Processed entry " << k + 1 << " / " << entries << '\n';
                }
            }

            std::cout << "[Group End] Finished processing group\n";
        }
        std::cout << "\n--- All Analysis Groups Finished ---\n";
        return true;
    }

    bool result() {
        if (m_reg.m_registry.empty()) {
            std::cout << "WARNING: Analysis registry is empty. Exiting result\n";
            return true;
        }

        for (auto const& [nav, analyses] : m_reg.m_registry) {
            for (const auto& analysis : analyses) {
                analysis->result();
            }
        }

        return true;
    }

    bool save() {
        if (m_reg.m_registry.empty()) {
            std::cout << "WARNING: Analysis registry is empty. Exiting result\n";
            return true;
        }

        for (auto const& [nav, analyses] : m_reg.m_registry) {
            for (const auto& analysis : analyses) {
                analysis->save();
            }
        }

        return true;
    }

private:

    analysis_registry& m_reg;

};

#endif // ANALYSIS_ANALYSIS_HPP_