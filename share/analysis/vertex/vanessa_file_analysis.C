#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TVector3.h>

struct prompt_event {
    int run_id;
    time_t sec;
    int nsec;
    double e;
};

bool operator<(const prompt_event& lhs, const prompt_event& rhs) {
    return TTimeStamp(lhs.sec, lhs.nsec) < TTimeStamp(rhs.sec, rhs.nsec);
}

struct ibd_debug {

    int run_id;

    time_t sec_p;
    int nsec_p;
    double e_p;
    double posx_p;
    double posy_p;
    double posz_p;

    time_t sec_d;
    int nsec_d;
    double e_d;
    double posx_d;
    double posy_d;
    double posz_d;

    bool prompt_cut;
    bool delayed_cut;
    bool fiducial_cut;
    bool chimney_cut;
    bool correlation_time_cut;
    bool correlation_space_cut;
    bool multiplicity_cut;
    bool neutron_cut;
    bool flasher_cut;

};

bool operator<(const ibd_debug& lhs, const ibd_debug& rhs) {
    return TTimeStamp(lhs.sec_p, lhs.nsec_p) < TTimeStamp(rhs.sec_p, rhs.nsec_p);
}

std::vector<double> generate_segment_boundaries(double start, double stop, int num_bins) {
    if (num_bins <= 0) return {};
    int num_points = num_bins + 1;
    double expected_width = (stop - start) / num_bins;
    
    std::vector<double> segment;
    segment.reserve(num_points);
    segment.push_back(start); 

    for (int i = 1; i < num_points; ++i) {
        double boundary = start + i * expected_width;
        
        if (i == num_points - 1) {
             segment.push_back(stop);
        } else {
             segment.push_back(boundary);
        }
    }
    return segment;
}

std::vector<double> create_custom_e_p_bins() {

    // double s1_start = 0.8;
    // double s2_start = 0.94;
    // double s3_start = 7.44;
    // double s4_start = 7.8;
    // double s5_start = 8.2;
    // double stop = 12.0;

    // int s1_bins = 1;
    // int s2_bins = 325;
    // int s3_bins = 9;
    // int s4_bins = 4;
    // int s5_bins = 1;
    // int tot_bins = s1_bins + s2_bins + s3_bins + s4_bins + s5_bins;

    double edges[] = {0.7, 1.0, 6.6, 7.4, 7.7, 8.1, 8.6, 9.4, 12.0};
    int    bins[]  = {  1,  56,   4,   1,   1,   1,   1,   1};
    
    std::vector<double> e_p_bins;

    for (std::size_t k = 0ul; k < 8ul; ++k) {
        double start = edges[k];
        double stop = edges[k + 1];
        int nbins = bins[k];
        std::vector<double> segment = generate_segment_boundaries(start, stop, nbins);
        if (k == 0ul) {
            e_p_bins.insert(e_p_bins.end(), segment.begin(), segment.end());
        } else {
            e_p_bins.insert(e_p_bins.end(), segment.begin() + 1, segment.end());
        }
    }
    
    return e_p_bins;
}

int vanessa_file_analysis(const char* filepath) {
    TFile* file = TFile::Open(filepath, "READ");
    if (!file) {
        std::cerr << "Error: Unable to open file " << filepath << '\n';
        return 1;
    }
    TTree* tree = file->Get<TTree>("ibd");
    if (!tree) {
        std::cerr << "Error: Cannot retrieve tree ibd in file " << filepath << '\n';
        return 1;
    }
    
    std::string* filename = nullptr;
    long long dt;
    long long dt_mu_time;
    double dist_track_time;
    int mu_type_time;
    int ntracks_time;
    int Nneu_time;
    char in_neu_veto_p;
    char in_neu_veto_d;
    long long dt_mu_dist;
    double dist_track_dist;
    int mu_type_dist;
    int ntracks_dist;
    int Nneu_dist;
    long long dt_mu_neu;
    double dist_track_neu;
    int mu_type_neu;
    int ntracks_neu;
    int Nneu_neu;
    double dR_p_neu;
    double dR_d_neu;
    std::vector<Long64_t>* mu_dt_all = nullptr;
    std::vector<int>* mu_type_all = nullptr;
    std::vector<int>* mu_Nneu_all = nullptr;
    std::vector<int>* mu_ntracks_all = nullptr;
    std::vector<double>* mu_dist_track_all = nullptr;
    std::vector<double>* mu_WP_PE_all = nullptr;
    std::vector<double>* mu_CD_PE_all = nullptr;
    double dR;
    double PE_p;
    double PE_d;
    double energy_p;
    double energy_d;
    double x_p;
    double y_p;
    double z_p;
    double x_d;
    double y_d;
    double z_d;
    double R_p;
    double R_d;
    double rho_p;
    double rho_d;
    long long t_p;
    long long t_d;
    int mult_after;
    int mult_between;
    int mult_before;
    int run;
    int trigid_p;
    int trigid_d;
    unsigned int oec_tag_p;
    unsigned int oec_tag_d;
    double energy_p_oec;
    double energy_d_oec;
    double x_p_oec;
    double y_p_oec;
    double z_p_oec;
    double x_d_oec;
    double y_d_oec;
    double z_d_oec;
    double dR_oec;
    char is_flasher_p;
    char is_lowEMuon_p;
    int nFired_p;
    float ntq_mean_p;
    float ntq_std_p;
    float hit_t_mean_p;
    float hit_t_std_p;
    float max_pmt_PE_p;
    float second_max_pmt_PE_p;
    float charge_ratio_p;
    float hit_q_mean_p;
    float hit_q_std_p;
    char is_flasher_d;
    char is_lowEMuon_d;
    int nFired_d;
    float ntq_mean_d;
    float ntq_std_d;
    float hit_t_mean_d;
    float hit_t_std_d;
    float max_pmt_PE_d;
    float second_max_pmt_PE_d;
    float charge_ratio_d;
    float hit_q_mean_d;
    float hit_q_std_d;

    tree->SetBranchAddress("filename", &filename);
    tree->SetBranchAddress("dt", &dt);
    tree->SetBranchAddress("dt_mu_time", &dt_mu_time);
    tree->SetBranchAddress("dist_track_time", &dist_track_time);
    tree->SetBranchAddress("mu_type_time", &mu_type_time);
    tree->SetBranchAddress("ntracks_time", &ntracks_time);
    tree->SetBranchAddress("Nneu_time", &Nneu_time);
    tree->SetBranchAddress("in_neu_veto_p", &in_neu_veto_p);
    tree->SetBranchAddress("in_neu_veto_d", &in_neu_veto_d);
    tree->SetBranchAddress("dt_mu_dist", &dt_mu_dist);
    tree->SetBranchAddress("dist_track_dist", &dist_track_dist);
    tree->SetBranchAddress("mu_type_dist", &mu_type_dist);
    tree->SetBranchAddress("ntracks_dist", &ntracks_dist);
    tree->SetBranchAddress("Nneu_dist", &Nneu_dist);
    tree->SetBranchAddress("dt_mu_neu", &dt_mu_neu);
    tree->SetBranchAddress("dist_track_neu", &dist_track_neu);
    tree->SetBranchAddress("mu_type_neu", &mu_type_neu);
    tree->SetBranchAddress("ntracks_neu", &ntracks_neu);
    tree->SetBranchAddress("Nneu_neu", &Nneu_neu);
    tree->SetBranchAddress("dR_p_neu", &dR_p_neu);
    tree->SetBranchAddress("dR_d_neu", &dR_d_neu);
    tree->SetBranchAddress("mu_dt_all", &mu_dt_all);
    tree->SetBranchAddress("mu_type_all", &mu_type_all);
    tree->SetBranchAddress("mu_Nneu_all", &mu_Nneu_all);
    tree->SetBranchAddress("mu_ntracks_all", &mu_ntracks_all);
    tree->SetBranchAddress("mu_dist_track_all", &mu_dist_track_all);
    tree->SetBranchAddress("mu_WP_PE_all", &mu_WP_PE_all);
    tree->SetBranchAddress("mu_CD_PE_all", &mu_CD_PE_all);
    tree->SetBranchAddress("dR", &dR);
    tree->SetBranchAddress("PE_p", &PE_p);
    tree->SetBranchAddress("PE_d", &PE_d);
    tree->SetBranchAddress("energy_p", &energy_p);
    tree->SetBranchAddress("energy_d", &energy_d);
    tree->SetBranchAddress("x_p", &x_p);
    tree->SetBranchAddress("y_p", &y_p);
    tree->SetBranchAddress("z_p", &z_p);
    tree->SetBranchAddress("x_d", &x_d);
    tree->SetBranchAddress("y_d", &y_d);
    tree->SetBranchAddress("z_d", &z_d);
    tree->SetBranchAddress("R_p", &R_p);
    tree->SetBranchAddress("R_d", &R_d);
    tree->SetBranchAddress("rho_p", &rho_p);
    tree->SetBranchAddress("rho_d", &rho_d);
    tree->SetBranchAddress("t_p", &t_p);
    tree->SetBranchAddress("t_d", &t_d);
    tree->SetBranchAddress("mult_after", &mult_after);
    tree->SetBranchAddress("mult_between", &mult_between);
    tree->SetBranchAddress("mult_before", &mult_before);
    tree->SetBranchAddress("run", &run);
    tree->SetBranchAddress("trigid_p", &trigid_p);
    tree->SetBranchAddress("trigid_d", &trigid_d);
    tree->SetBranchAddress("oec_tag_p", &oec_tag_p);
    tree->SetBranchAddress("oec_tag_d", &oec_tag_d);
    tree->SetBranchAddress("energy_p_oec", &energy_p_oec);
    tree->SetBranchAddress("energy_d_oec", &energy_d_oec);
    tree->SetBranchAddress("x_p_oec", &x_p_oec);
    tree->SetBranchAddress("y_p_oec", &y_p_oec);
    tree->SetBranchAddress("z_p_oec", &z_p_oec);
    tree->SetBranchAddress("x_d_oec", &x_d_oec);
    tree->SetBranchAddress("y_d_oec", &y_d_oec);
    tree->SetBranchAddress("z_d_oec", &z_d_oec);
    tree->SetBranchAddress("dR_oec", &dR_oec);
    tree->SetBranchAddress("is_flasher_p", &is_flasher_p);
    tree->SetBranchAddress("is_lowEMuon_p", &is_lowEMuon_p);
    tree->SetBranchAddress("nFired_p", &nFired_p);
    tree->SetBranchAddress("ntq_mean_p", &ntq_mean_p);
    tree->SetBranchAddress("ntq_std_p", &ntq_std_p);
    tree->SetBranchAddress("hit_t_mean_p", &hit_t_mean_p);
    tree->SetBranchAddress("hit_t_std_p", &hit_t_std_p);
    tree->SetBranchAddress("max_pmt_PE_p", &max_pmt_PE_p);
    tree->SetBranchAddress("second_max_pmt_PE_p", &second_max_pmt_PE_p);
    tree->SetBranchAddress("charge_ratio_p", &charge_ratio_p);
    tree->SetBranchAddress("hit_q_mean_p", &hit_q_mean_p);
    tree->SetBranchAddress("hit_q_std_p", &hit_q_std_p);
    tree->SetBranchAddress("is_flasher_d", &is_flasher_d);
    tree->SetBranchAddress("is_lowEMuon_d", &is_lowEMuon_d);
    tree->SetBranchAddress("nFired_d", &nFired_d);
    tree->SetBranchAddress("ntq_mean_d", &ntq_mean_d);
    tree->SetBranchAddress("ntq_std_d", &ntq_std_d);
    tree->SetBranchAddress("hit_t_mean_d", &hit_t_mean_d);
    tree->SetBranchAddress("hit_t_std_d", &hit_t_std_d);
    tree->SetBranchAddress("max_pmt_PE_d", &max_pmt_PE_d);
    tree->SetBranchAddress("second_max_pmt_PE_d", &second_max_pmt_PE_d);
    tree->SetBranchAddress("charge_ratio_d", &charge_ratio_d);
    tree->SetBranchAddress("hit_q_mean_d", &hit_q_mean_d);
    tree->SetBranchAddress("hit_q_std_d", &hit_q_std_d);

    std::vector<double> e_p_bins = create_custom_e_p_bins();
    TH1D* h_e_p = new TH1D("h_e_p", "Prompt energy;E_{p} (MeV);Entries;", e_p_bins.size() - 1, e_p_bins.data());

    int run_min = 999999, run_max = 0;

    // std::set<prompt_event> events;
    std::set<ibd_debug> events;

    for (long k = 0l; k < tree->GetEntries(); ++k) {
        tree->GetEntry(k);

        if (run < 9789 || 9822 < run) continue;

        if (energy_p < 0.7 || 12.0 < energy_p) continue;
        if (energy_d < 2.0 || 2.5 < energy_d) continue;
        ibd_debug ibd{
            run,
            static_cast<time_t>(t_p / 1000000000),
            static_cast<int>(t_p % 1000000000),
            energy_p,
            x_p,
            y_p,
            z_p,
            static_cast<time_t>(t_d / 1000000000),
            static_cast<int>(t_d % 1000000000),
            energy_d,
            x_d,
            y_d,
            z_d,
            0.7 <= energy_p && energy_p <= 12.0,
            2.0 <= energy_d && energy_d <= 2.5,
            R_p <= 16.5,
            std::abs(z_p) <= 15.5 || rho_p >= 2.0,
            5000 <= dt && dt <= 1000000,
            dR <= 1.5,
            !mult_before && !mult_between && !mult_after,
            !in_neu_veto_p && !in_neu_veto_d,
            !is_flasher_p && !is_flasher_d
        };
        events.insert(ibd);
        
        // if (energy_p < 0.7 || 12.0 < energy_p) continue;
        // if ((energy_d < 2.0 || 2.5 < energy_d) && (energy_d < 4.5 || 5.5 < energy_d)) continue;
        // if (16.5 < R_p) continue;
        // if (std::abs(z_p) > 15.5 && rho_p < 2.0) continue;
        // if (dt < 5000 || 1000000 < dt) continue;
        // if (1.5 < dR) continue;
        // if (mult_before || mult_between || mult_after) continue;
        // if (in_neu_veto_p || in_neu_veto_d) continue;
        // if (is_flasher_p || is_flasher_d) continue;
        // events.insert(
        //     prompt_event{run, static_cast<time_t>(t_p / 1000000000), static_cast<int>(t_p % 1000000000), energy_p}
        // );
        // h_e_p->Fill(energy_p);
    }

    for (const ibd_debug& ibd : events) {
        std::cout << ibd.run_id << ' ' 
                  << ibd.sec_p << ' ' << ibd.nsec_p << ' ' << ibd.e_p << ' ' << ibd.posx_p << ' ' << ibd.posy_p << ' ' << ibd.posz_p << ' ' 
                  << ibd.sec_d << ' ' << ibd.nsec_d << ' ' << ibd.e_d << ' ' << ibd.posx_d << ' ' << ibd.posy_d << ' ' << ibd.posz_d << ' '
                  << ibd.prompt_cut << ' ' << ibd.delayed_cut << ' ' << ibd.fiducial_cut << ' ' << ibd.chimney_cut << ' ' 
                  << ibd.correlation_time_cut << ' ' << ibd.correlation_space_cut << ' ' 
                  << ibd.multiplicity_cut << ' ' << ibd.neutron_cut << ' ' << ibd.flasher_cut << '\n';
    }

    // for (const prompt_event& prompt : events) {
    //     std::cout << prompt.run_id << ' ' << prompt.sec << ' ' << prompt.nsec << ' ' << prompt.e << '\n';
    // }

    TCanvas* c_e_p = new TCanvas("c_e_p", "c_e_p", 1000, 1000);
    c_e_p->cd();

    h_e_p->GetXaxis()->CenterTitle(true);
    h_e_p->GetYaxis()->CenterTitle(true);
    h_e_p->GetXaxis()->SetTitleOffset(1.25);
    h_e_p->SetLineWidth(2);
    h_e_p->SetLineColor(kBlue);
    h_e_p->Draw("HIST");

    c_e_p->SetTickx();
    c_e_p->SetTicky();
    c_e_p->Update();

    return 0;
}