#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>

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
    
    std::string filename;
    double dt;
    double dt_mu_time;
    double dist_track_time;
    int mu_type_time;
    int ntracks_time;
    int Nneu_time;
    bool in_neu_veto_p;
    bool in_neu_veto_d;
    double dt_mu_dist;
    double dist_track_dist;
    int mu_type_dist;
    int ntracks_dist;
    int Nneu_dist;
    double dt_mu_neu;
    double dist_track_neu;
    int mu_type_neu;
    int ntracks_neu;
    int Nneu_neu;
    double dR_p_neu;
    double dR_d_neu;
    std::vector<Long64_t> mu_dt_all;
    std::vector<int> mu_type_all;
    std::vector<int> mu_Nneu_all;
    std::vector<int> mu_ntracks_all;
    std::vector<double> mu_dist_track_all;
    std::vector<double> mu_WP_PE_all;
    std::vector<double> mu_CD_PE_all;
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
    double t_p;
    double t_d;
    int mult_after;
    int mult_between;
    int mult_before;
    int run;
    int trigid_p;
    int trigid_d;
    int oec_tag_p;
    int oec_tag_d;
    double energy_p_oec;
    double energy_d_oec;
    double x_p_oec;
    double y_p_oec;
    double z_p_oec;
    double x_d_oec;
    double y_d_oec;
    double z_d_oec;
    double dR_oec;
    bool is_flasher_p;
    bool is_lowEMuon_p;
    int nFired_p;
    double ntq_mean_p;
    double ntq_std_p;
    double hit_t_mean_p;
    double hit_t_std_p;
    double max_pmt_PE_p;
    double second_max_pmt_PE_p;
    double charge_ratio_p;
    double hit_q_mean_p;
    double hit_q_std_p;
    bool is_flasher_d;
    bool is_lowEMuon_d;
    int nFired_d;
    double ntq_mean_d;
    double ntq_std_d;
    double hit_t_mean_d;
    double hit_t_std_d;
    double max_pmt_PE_d;
    double second_max_pmt_PE_d;
    double charge_ratio_d;
    double hit_q_mean_d;
    double hit_q_std_d;

    return 0;
}