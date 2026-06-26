from datetime import datetime, timestamp

import numpy as np

def apply_veto(ts_target : np.ndarray, ts_all : np.ndarray, veto_window : timestamp = timestamp(0, 10000000)) -> np.ndarray:
    mask_keep = np.ones(len(ts_target), dtype=bool)
    veto_window_full = np.full(len(ts_all), veto_window)
    for i, ts in enumerate(ts_target):
        dt = ts_all - ts
        if np.any((veto_window_full < dt) & (dt < veto_window_full) & (dt != timestamp())):
            mask_keep[i] = False
    return ts_target[mask_keep]

def calculate_wp_tagging_efficiency(ts_all : np.ndarray, ts_cd_only : np.ndarray, rate_total_cd : float, rate_total_cd_err : float, nb_cd_wp_high : int):
    # ⬅️ = avant processing
    # ➡️ = après processing
    # ⬆️ = paramètre

    # Calculation

    # rate_total ➡️✅ = rate_cd_only ➡️✅ + rate_cd_wp ➡️✅
    # veto_window ⬆️✅

    # nb_cd_wp_high ⬅️✅

    veto_window = timestamp(0, 10000000)

    prob_veto = 1.0 - np.exp(-rate_total_cd * veto_window.to_sec())

    nb_cd_only_veto = len(apply_veto(ts_cd_only, ts_all, veto_window))
    nb_cd_only_corr = nb_cd_only_veto * (1.0 - prob_veto)
    wp_tagging_efficiency = 1.0 - (nb_cd_only_corr / nb_cd_wp_high)

    # Error propogation ✅

    nb_cd_only_veto_err = np.sqrt(nb_cd_only_veto)
    prob_veto_err = veto_window.to_sec() * np.exp(-rate_total_cd * veto_window.to_sec()) * rate_total_cd_err
    nb_cd_only_corr_err = np.sqrt(
        ((1 - prob_veto) * nb_cd_only_veto_err)**2 + 
        (nb_cd_only_veto * prob_veto_err)**2
    )
    nb_cd_wp_high_err = np.sqrt(nb_cd_wp_high)
    wp_tagging_efficiency_err = np.sqrt(
        (nb_cd_only_corr_err / nb_cd_wp_high)**2 +
        (nb_cd_only_corr * nb_cd_wp_high_err / nb_cd_wp_high**2)**2
    )

    return wp_tagging_efficiency, wp_tagging_efficiency_err

def calculate_muon_rate(run : int, plot=False):

    data = {}
    mask_cd_only = np.logical_and(data["totq_cd"] > 0, data["totq_wp"] == 0)
    mask_wp_only = np.logical_and(data["totq_cd"] == 0, data["totq_wp"] > 0)

    data_cd_only = {key: val[mask_cd_only] for key, val in data.items()}
    data_wp_only = {key: val[mask_wp_only] for key, val in data.items()}

    ts_cd_only = np.array([timestamp(sec, nsec) for sec, nsec in zip(data_cd_only["sec"], data_cd_only["nsec"])])
    ts_wp_only = np.array([timestamp(sec, nsec) for sec, nsec in zip(data_wp_only["sec"], data_wp_only["nsec"])])

    ts_all = np.concatenate([ts_cd_wp, ts_cd_only, ts_wp_only])
    wp_tagging_efficiency, wp_tagging_efficiency_err = calculate_wp_tagging_efficiency(ts_all, ts_cd_only)