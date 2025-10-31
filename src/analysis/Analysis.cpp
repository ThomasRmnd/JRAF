#include "analysis/Analysis.hpp"

#include <iostream>

#include "event/Event.hpp"
#include "event/EventCache.hpp"

Analysis::Analysis(const std::string& name, const std::string& method) :
    m_name{name},
    m_method_sel{method}
{}

bool Analysis::initialize() {
    m_tree = new TTree(m_name.c_str(), m_name.c_str());
    if (!m_tree) {
        std::cout << "[ERROR] Could not create the TTree for Analysis " << m_name << '\n';
        return false;
    }

    m_tree->Branch("posx_p", &posx_p);
    m_tree->Branch("posy_p", &posy_p);
    m_tree->Branch("posz_p", &posz_p);
    m_tree->Branch("e_p", &e_p);
    m_tree->Branch("totq_p", &totq_p);
    m_tree->Branch("sec_p", &sec_p);
    m_tree->Branch("nsec_p", &nsec_p);

    m_tree->Branch("posx_d", &posx_d);
    m_tree->Branch("posy_d", &posy_d);
    m_tree->Branch("posz_d", &posz_d);
    m_tree->Branch("e_d", &e_d);
    m_tree->Branch("totq_d", &totq_d);
    m_tree->Branch("sec_d", &sec_d);
    m_tree->Branch("nsec_d", &nsec_d);
    return true;
}

bool Analysis::write() {
    if (m_tree) {
        m_tree->Write();
    }
    return true;
}

void Analysis::extractEvent(JM::NavBuffer* buf, std::vector<std::vector<track>>& tracks, std::vector<vertex>& cur_vertices, std::vector<vertex>& bef_vertices, std::vector<vertex>& aft_vertices) {
    tracks.reserve(buf->size());
    bef_vertices.reserve(buf->size() / 2);
    aft_vertices.reserve(buf->size() / 2);

    for (JM::NavBuffer::Iterator it = buf->begin(); it != buf->end(); ++it) {
        JM::EvtNavigator* nav = it->get();
        if (!nav) continue;

        std::shared_ptr<Event> evt_ptr = EventCache::load(nav);
        if (!evt_ptr) continue;

        const Event& evt = *evt_ptr;

        tracks.push_back(evt.tracks);

        if (it < buf->current()) {
            for (const auto& v : evt.vertices) {
                if (m_method_sel.isIn(v)) {
                    bef_vertices.push_back(v);
                }
            }
        } else if (buf->current() < it) {
            for (const auto& v : evt.vertices) {
                if (m_method_sel.isIn(v)) {
                    aft_vertices.push_back(v);
                }
            }
        } else {
            for (const auto& v : evt.vertices) {
                if (m_method_sel.isIn(v)) {
                    cur_vertices.push_back(v);
                }
            }
        }
    }
}