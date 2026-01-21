if (bool_perc == "TRUE") {

        // TODO: add the 68.3% percentile for the cylindrical volume ratio

        std::vector<std::vector<double>> vecs = {v_angle, v_dist_mid_point, v_min_dist_from_init_point, v_min_dist_from_end_point, v_dist_bottom_wp};
        double max_len = 18000.0 * 18000.0;
        const std::size_t nb_bins = 9;
        std::vector<TH1D*> graphs_68perc{vecs.size(), nullptr};
        std::vector<TCanvas*> canvases_68perc{vecs.size(), nullptr};

        std::vector<double> dist_center;
        for (int i = 0; i < tree->GetEntries(); ++i) {
            tree->GetEntry(i);
            dist_center.push_back(y_axis);
        }

        for (std::size_t k = 0; k < vecs.size(); ++k) {
            graphs_68perc[k] = new TH1D(("graphs_68perc__" + std::to_string(k)).c_str(), "", nb_bins, 0.0, max_len);
            graphs_68perc[k]->SetStats(kFALSE);
            graphs_68perc[k]->SetMarkerStyle(kFullCircle);
            graphs_68perc[k]->SetMarkerSize(2);
            graphs_68perc[k]->SetMarkerColor(k + 1);
            graphs_68perc[k]->SetLineWidth(2);
            graphs_68perc[k]->SetLineColor(k + 1);
            graphs_68perc[k]->GetXaxis()->SetNdivisions(nb_bins, kFALSE);
            std::array<std::vector<double>, nb_bins> vecs_68perc;
            for (std::size_t i = 0; i < dist_center.size(); ++i) {
                std::size_t j = std::floor((dist_center[i] * dist_center[i]) * nb_bins / max_len);
                vecs_68perc[j].push_back(vecs[k][i]);
            }
            for (std::size_t i = 0; i < nb_bins; ++i) {
                std::cout << "[" << std::sqrt((max_len * i) / nb_bins) << "^2, " << std::sqrt((max_len * (i + 1)) / nb_bins) << "^2[: " << vecs_68perc[i].size() << std::endl;
                if (vecs_68perc[i].size() == 0) {
                    graphs_68perc[k]->SetBinContent(i + 1, 0.0);
                    continue;
                }
                std::nth_element(vecs_68perc[i].begin(), vecs_68perc[i].begin() + vecs_68perc[i].size() * 683 / 1000, vecs_68perc[i].end());
                graphs_68perc[k]->SetBinContent(i + 1, vecs_68perc[i][vecs_68perc[i].size() * 683 / 1000]);
            }
            for (std::size_t i = 0; i < nb_bins; ++i) {
                graphs_68perc[k]->SetBinError(i + 1, 0.00001);
                std::stringstream ss;
                ss << std::fixed << std::setprecision(1) << std::sqrt((max_len * i) / nb_bins) / 1000.0;
                graphs_68perc[k]->GetXaxis()->ChangeLabel(i + 1, -1.0, -1.0, -1, -1, -1, (ss.str() + "^{2}").c_str());
            }
            std::stringstream ss;
            ss << std::fixed << std::setprecision(1) << std::sqrt(max_len) / 1000.0;
            graphs_68perc[k]->GetXaxis()->ChangeLabel(nb_bins + 1, -1.0, -1.0, -1, -1, -1, (ss.str() + "^{2}").c_str());
            canvases_68perc[k] = new TCanvas(("c_68perc_" + std::to_string(k)).c_str(), ("c" + std::to_string(k + 6)).c_str(), 1200, 1200);
            canvases_68perc[k]->cd();
            canvases_68perc[k]->SetGrid();
            graphs_68perc[k]->Draw("E P SAMES");
            canvases_68perc[k]->Update();
        }

        for (std::vector<double>& vec : vecs) {
            std::nth_element(vec.begin(), vec.begin() + vec.size() * 683 / 1000, vec.end());
            std::cout << "68.3% quantile: " << vec[vec.size() * 683 / 1000] << std::endl;
            std::nth_element(vec.begin(), vec.begin() + vec.size() * 955 / 1000, vec.end());
            std::cout << "95.5% quantile: " << vec[vec.size() * 955 / 1000] << std::endl;
        }

    }