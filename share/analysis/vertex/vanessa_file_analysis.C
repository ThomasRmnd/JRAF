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
    std::vector<double> dR;
    std::vector<double> PE_p;
    std::vector<double> PE_d;
    std::vector<double> energy_p;
    std::vector<double> energy_d;
    std::vector<double> x_p;
    std::vector<double> y_p;
    std::vector<double> z_p;
    std::vector<double> x_d;
    std::vector<double> y_d;
    std::vector<double> z_d;
    std::vector<double> R_p;
    std::vector<double> R_d;
    std::vector<double> rho_p;
    std::vector<double> rho_d;


*Br   39 :R_p       : R_p/D                                                  *
*Entries :    14857 : Total  Size=     535681 bytes  File Size  =     441776 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   40 :R_d       : R_d/D                                                  *
*Entries :    14857 : Total  Size=     535681 bytes  File Size  =     441776 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   41 :rho_p     : rho_p/D                                                *
*Entries :    14857 : Total  Size=     545049 bytes  File Size  =     451136 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   42 :rho_d     : rho_d/D                                                *
*Entries :    14857 : Total  Size=     545049 bytes  File Size  =     451136 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   43 :t_p       : t_p/L                                                  *
*Entries :    14857 : Total  Size=     535681 bytes  File Size  =     441776 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   44 :t_d       : t_d/L                                                  *
*Entries :    14857 : Total  Size=     535681 bytes  File Size  =     441773 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   45 :mult_after : mult_after/I                                          *
*Entries :    14857 : Total  Size=     509033 bytes  File Size  =     411999 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.01     *
*Br   46 :mult_between : mult_between/I                                      *
*Entries :    14857 : Total  Size=     518401 bytes  File Size  =     421283 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.01     *
*Br   47 :mult_before : mult_before/I                                        *
*Entries :    14857 : Total  Size=     513717 bytes  File Size  =     416696 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.01     *
*Br   48 :run       : run/I                                                  *
*Entries :    14857 : Total  Size=     476245 bytes  File Size  =     380766 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   49 :trigid_p  : trigid_p/I                                             *
*Entries :    14857 : Total  Size=     499665 bytes  File Size  =     405748 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   50 :trigid_d  : trigid_d/I                                             *
*Entries :    14857 : Total  Size=     499665 bytes  File Size  =     405748 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   51 :oec_tag_p : oec_tag_p/i                                            *
*Entries :    14857 : Total  Size=     504349 bytes  File Size  =     410404 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   52 :oec_tag_d : oec_tag_d/i                                            *
*Entries :    14857 : Total  Size=     504349 bytes  File Size  =     410362 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   53 :energy_p_oec : energy_p_oec/D                                      *
*Entries :    14857 : Total  Size=     577837 bytes  File Size  =     483869 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   54 :energy_d_oec : energy_d_oec/D                                      *
*Entries :    14857 : Total  Size=     577837 bytes  File Size  =     483745 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   55 :x_p_oec   : x_p_oec/D                                              *
*Entries :    14857 : Total  Size=     554417 bytes  File Size  =     460496 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   56 :y_p_oec   : y_p_oec/D                                              *
*Entries :    14857 : Total  Size=     554417 bytes  File Size  =     460496 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   57 :z_p_oec   : z_p_oec/D                                              *
*Entries :    14857 : Total  Size=     554417 bytes  File Size  =     460496 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   58 :x_d_oec   : x_d_oec/D                                              *
*Entries :    14857 : Total  Size=     554417 bytes  File Size  =     460496 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   59 :y_d_oec   : y_d_oec/D                                              *
*Entries :    14857 : Total  Size=     554417 bytes  File Size  =     460496 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   60 :z_d_oec   : z_d_oec/D                                              *
*Entries :    14857 : Total  Size=     554417 bytes  File Size  =     460496 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   61 :dR_oec    : dR_oec/D                                               *
*Entries :    14857 : Total  Size=     549733 bytes  File Size  =     455816 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   62 :is_flasher_p : is_flasher_p/B                                      *
*Entries :    14857 : Total  Size=     473824 bytes  File Size  =     379897 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   63 :is_lowEMuon_p : is_lowEMuon_p/B                                    *
*Entries :    14857 : Total  Size=     478508 bytes  File Size  =     384577 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   64 :nFired_p  : nFired_p/I                                             *
*Entries :    14857 : Total  Size=     499665 bytes  File Size  =     405748 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   65 :ntq_mean_p : ntq_mean_p/F                                          *
*Entries :    14857 : Total  Size=     509033 bytes  File Size  =     415108 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   66 :ntq_std_p : ntq_std_p/F                                            *
*Entries :    14857 : Total  Size=     504349 bytes  File Size  =     410428 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   67 :hit_t_mean_p : hit_t_mean_p/F                                      *
*Entries :    14857 : Total  Size=     518401 bytes  File Size  =     424468 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   68 :hit_t_std_p : hit_t_std_p/F                                        *
*Entries :    14857 : Total  Size=     513717 bytes  File Size  =     419788 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   69 :max_pmt_PE_p : max_pmt_PE_p/F                                      *
*Entries :    14857 : Total  Size=     518401 bytes  File Size  =     424468 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   70 :second_max_pmt_PE_p : second_max_pmt_PE_p/F                        *
*Entries :    14857 : Total  Size=     551189 bytes  File Size  =     457228 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   71 :charge_ratio_p : charge_ratio_p/F                                  *
*Entries :    14857 : Total  Size=     527769 bytes  File Size  =     433828 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   72 :hit_q_mean_p : hit_q_mean_p/F                                      *
*Entries :    14857 : Total  Size=     518401 bytes  File Size  =     424468 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   73 :hit_q_std_p : hit_q_std_p/F                                        *
*Entries :    14857 : Total  Size=     513717 bytes  File Size  =     419788 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   74 :is_flasher_d : is_flasher_d/B                                      *
*Entries :    14857 : Total  Size=     473824 bytes  File Size  =     379897 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   75 :is_lowEMuon_d : is_lowEMuon_d/B                                    *
*Entries :    14857 : Total  Size=     478508 bytes  File Size  =     384577 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   76 :nFired_d  : nFired_d/I                                             *
*Entries :    14857 : Total  Size=     499665 bytes  File Size  =     405748 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   77 :ntq_mean_d : ntq_mean_d/F                                          *
*Entries :    14857 : Total  Size=     509033 bytes  File Size  =     415108 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   78 :ntq_std_d : ntq_std_d/F                                            *
*Entries :    14857 : Total  Size=     504349 bytes  File Size  =     410428 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   79 :hit_t_mean_d : hit_t_mean_d/F                                      *
*Entries :    14857 : Total  Size=     518401 bytes  File Size  =     424468 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   80 :hit_t_std_d : hit_t_std_d/F                                        *
*Entries :    14857 : Total  Size=     513717 bytes  File Size  =     419788 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   81 :max_pmt_PE_d : max_pmt_PE_d/F                                      *
*Entries :    14857 : Total  Size=     518401 bytes  File Size  =     424468 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   82 :second_max_pmt_PE_d : second_max_pmt_PE_d/F                        *
*Entries :    14857 : Total  Size=     551189 bytes  File Size  =     457228 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   83 :charge_ratio_d : charge_ratio_d/F                                  *
*Entries :    14857 : Total  Size=     527769 bytes  File Size  =     433828 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   84 :hit_q_mean_d : hit_q_mean_d/F                                      *
*Entries :    14857 : Total  Size=     518401 bytes  File Size  =     424468 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *
*Br   85 :hit_q_std_d : hit_q_std_d/F                                        *
*Entries :    14857 : Total  Size=     513717 bytes  File Size  =     419788 *
*Baskets :     4680 : Basket Size=      32000 bytes  Compression=   1.00     *


    return 0;
}