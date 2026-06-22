import argparse
import os

from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import uproot

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ibd-analysis", type=str, nargs="+", help="IBD analysis filepath")
    parser.add_argument("--cosmo-shape-analysis", type=str, nargs="+", help="Cosmo shape analysis filepath")
    parser.add_argument("--run-info", type=str, help="Run info filepath")
    return parser.parse_args()

def cosmo_shape_analysis_plot(filepath: str, **meta):
    file = uproot.open(filepath)
    tree_bkg = file["background_events"]
    tree_sig = file["signal_events"]

    branches = [
        "run_id",
        "posx_p", "posy_p", "posz_p", "sec_p", "nsec_p", "e_p", "dlat_mu2p", "dt_mu2p",
        "posx_d", "posy_d", "posz_d", "sec_d", "nsec_d", "e_d", "dlat_mu2d", "dt_mu2d"
    ]

    data_bkg = tree_bkg.arrays(branches, library="np")
    data_sig = tree_sig.arrays(branches, library="np")

    # mask_sig_reprodc = np.logical_and(9789 <= data_sig["run_id"], data_sig["run_id"] <= 11039)
    # mask_sig_reprodd = np.logical_and(11049 <= data_sig["run_id"], data_sig["run_id"] <= 12135)
    # data_sig = {key: value[mask_sig_reprodd] for key, value in data_sig.items()}
    # mask_bkg_reprodc = np.logical_and(9789 <= data_bkg["run_id"], data_bkg["run_id"] <= 11039)
    # mask_sig_reprodd = np.logical_and(11049 <= data_bkg["run_id"], data_bkg["run_id"] <= 12135)
    # data_bkg = {key: value[mask_sig_reprodd] for key, value in data_bkg.items()}

    ts_p_bkg = np.array([timestamp(sec, nsec) for sec, nsec in zip(data_bkg["sec_p"], data_bkg["nsec_p"])])
    ts_d_bkg = np.array([timestamp(sec, nsec) for sec, nsec in zip(data_bkg["sec_d"], data_bkg["nsec_d"])])
    ts_diff_bkg = np.array([(t_d - t_p).to_sec() * 1e3 for t_p, t_d in zip(ts_p_bkg, ts_d_bkg)]) # in ms

    pos_p_bkg = np.column_stack((data_bkg["posx_p"], data_bkg["posy_p"], data_bkg["posz_p"]))
    pos_d_bkg = np.column_stack((data_bkg["posx_d"], data_bkg["posy_d"], data_bkg["posz_d"]))
    distance_bkg = np.linalg.norm(pos_p_bkg - pos_d_bkg, axis=1) / 1000.0
    rho_p_bkg = np.linalg.norm(pos_p_bkg[:, :2] / 1000.0, axis=1) ** 2
    z_p_bkg = pos_p_bkg[:, 2] / 1000.0
    rho_d_bkg = np.linalg.norm(pos_d_bkg[:, :2] / 1000.0, axis=1) ** 2
    z_d_bkg = pos_d_bkg[:, 2] / 1000.0

    ts_p_sig = np.array([timestamp(sec, nsec) for sec, nsec in zip(data_sig["sec_p"], data_sig["nsec_p"])])
    ts_d_sig = np.array([timestamp(sec, nsec) for sec, nsec in zip(data_sig["sec_d"], data_sig["nsec_d"])])
    ts_diff_sig = np.array([(t_d - t_p).to_sec() * 1e3 for t_p, t_d in zip(ts_p_sig, ts_d_sig)]) # in ms

    pos_p_sig = np.column_stack((data_sig["posx_p"], data_sig["posy_p"], data_sig["posz_p"]))
    pos_d_sig = np.column_stack((data_sig["posx_d"], data_sig["posy_d"], data_sig["posz_d"]))
    distance_sig = np.linalg.norm(pos_p_sig - pos_d_sig, axis=1) / 1000.0
    rho_p_sig = np.linalg.norm(pos_p_sig[:, :2] / 1000.0, axis=1) ** 2
    z_p_sig = pos_p_sig[:, 2] / 1000.0
    rho_d_sig = np.linalg.norm(pos_d_sig[:, :2] / 1000.0, axis=1) ** 2
    z_d_sig = pos_d_sig[:, 2] / 1000.0

    print(f"Loaded {len(data_bkg['posx_p'])} background events in {filepath}")
    print(f"Loaded {len(data_sig['posx_p'])} signal events in {filepath}")

    mask_reprod25d_sig = data_sig["run_id"] >= 11049
    mask_reprod25d_bkg = data_bkg["run_id"] >= 11049
    data_sig["e_p"][mask_reprod25d_sig] *= 1.00950656406
    data_sig["e_d"][mask_reprod25d_sig] *= 1.00950656406
    data_bkg["e_p"][mask_reprod25d_bkg] *= 1.00950656406
    data_bkg["e_d"][mask_reprod25d_bkg] *= 1.00950656406

    e_p_plotter = PromptEnergyPlotter(binmode="nmo")
    e_p_plotter.add(data_sig["e_p"], linecolor="#648fff", fillcolor="#eff3ff", label="Cosmogenic enriched")
    e_p_plotter.add(data_bkg["e_p"], linecolor="#ff6464", fillcolor="#ffefef", label="Cosmogenic depleted") 
    e_p_plotter.plot(**meta) 
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_e_p_nmo.pdf')}") 
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_e_p_nmo.png')}")

    e_p_uncertainty_plotter = PromptEnergyUncertaintyPlotter(binmode="normal")
    e_p_uncertainty_plotter.add(data_sig["e_p"], linecolor="#648fff", fillcolor="#eff3ff", label="Cosmogenic enriched")
    e_p_uncertainty_plotter.add(data_bkg["e_p"], linecolor="#ff6464", fillcolor="#ffefef", label="Cosmogenic depleted")
    e_p_uncertainty_plotter.plot(**meta)

    e_p_plotter_normal = PromptEnergyPlotter(binmode="normal", force_not_plot_main=True, mc_cosmo_spectrum="mc/mc_cosmogenics.root")
    e_p_plotter_normal.add(data_sig["e_p"], linecolor="#648fff", label="Cosmogenic enriched")
    e_p_plotter_normal.add(data_bkg["e_p"], linecolor="#ff6464", label="Cosmogenic depleted")
    e_p_plotter_normal.plot(**meta)

    e_p_plotter_normal = PromptEnergyPlotter(binmode="normal")
    e_p_plotter_normal.add(data_sig["e_p"], linecolor="#648fff", label="Cosmogenic enriched")
    e_p_plotter_normal.add(data_bkg["e_p"], linecolor="#ff6464", label="Cosmogenic depleted")
    e_p_plotter_normal.plot(**meta)
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_e_p_normal.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_e_p_normal.png')}")

    e_d_plotter = DelayedEnergyPlotter()
    e_d_plotter.add(data_sig["e_d"], linecolor="#648fff", label="Cosmogenic enriched")
    e_d_plotter.add(data_bkg["e_d"], linecolor="#ff6464", label="Cosmogenic depleted")
    e_d_plotter.plot()
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_e_d.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_e_d.png')}")

    ts_diff_plotter = PromptDelayedTimePlotter()
    ts_diff_plotter.add(ts_diff_sig, linecolor="#648fff", label="Cosmogenic enriched")
    ts_diff_plotter.add(ts_diff_bkg, linecolor="#ff6464", label="Cosmogenic depleted")
    ts_diff_plotter.plot()
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_dt_p_d.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_dt_p_d.png')}")

    distance_plotter = PromptDelayedDistancePlotter()
    distance_plotter.add(distance_sig, linecolor="#648fff", label="Cosmogenic enriched")
    distance_plotter.add(distance_bkg, linecolor="#ff6464", label="Cosmogenic depleted")
    distance_plotter.plot()
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_dr_p_d.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_dr_p_d.png')}")

    spatial_plotter = SpatialDistributionPlotter()
    spatial_plotter.plot(rho_p_bkg, z_p_bkg, r"$\rho_{p}$ (m$^2$)", r"$z_{p}$ (m)")
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_rho_z_p_bkg.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_rho_z_p_bkg.png')}")
    spatial_plotter.plot(rho_d_bkg, z_d_bkg, r"$\rho_{d}$ (m$^2$)", r"$z_{d}$ (m)")
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_rho_z_d_bkg.pdf')}") 
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_rho_z_d_bkg.png')}")
    spatial_plotter.plot(rho_p_sig, z_p_sig, r"$\rho_{p}$ (m$^2$)", r"$z_{p}$ (m)")
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_rho_z_p_sig.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_rho_z_p_sig.png')}")
    spatial_plotter.plot(rho_d_sig, z_d_sig, r"$\rho_{d}$ (m$^2$)", r"$z_{d}$ (m)")
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_rho_z_d_sig.pdf')}") 
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_rho_z_d_sig.png')}")

    xbins_bkg = np.linspace(-1.2, 0.0, 51)
    xbins_sig = np.linspace(0.0, 1.2, 51)
    if "_3m_" in filepath:
        ybins = np.linspace(0.0, 3.0, 51)
    elif "_4m_" in filepath:
        ybins = np.linspace(0.0, 4.0, 51)
    elif "_1_5m_" in filepath:
        ybins = np.linspace(0.0, 1.5, 51)

    if "_1_2s_" in filepath:
        xbins_bkg = np.linspace(-1.2, 0.0, 51)
        xbins_sig = np.linspace(0.0, 1.2, 51)
    elif "_2s_" in filepath:
        xbins_bkg = np.linspace(-2.0, 0.0, 51)
        xbins_sig = np.linspace(0.0, 2.0, 51)

    muon_veto_plotter = MuonVetoDistributionPlotter()
    muon_veto_plotter.plot(data_sig["dt_mu2p"], data_sig["dlat_mu2p"], xbins_sig, ybins, r"$\Delta t_{\mu-p}$ (s)", r"$d_{\mu-p}$ (m)", is_signal_region=True)
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_dt_dlat_p_sig.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_dt_dlat_p_sig.png')}")
    muon_veto_plotter.plot(data_bkg["dt_mu2p"], data_bkg["dlat_mu2p"], xbins_bkg, ybins, r"$\Delta t_{\mu-p}$ (s)", r"$d_{\mu-p}$ (m)", is_signal_region=False)
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_dt_dlat_p_bkg.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_dt_dlat_p_bkg.png')}")

    plt.show()

def ibd_analysis_plot(filepath: str, **meta):
    file = uproot.open(filepath)
    tree = file["events"]

    branches = [
        "run_id",
        "posx_p", "posy_p", "posz_p", "sec_p", "nsec_p", "e_p",
        "posx_d", "posy_d", "posz_d", "sec_d", "nsec_d", "e_d",
        "dt_last_mu", "dt_mu2p", "dlat_mu2p"
    ]
    data = tree.arrays(branches, library="np")

    ts_p = np.array([timestamp(sec, nsec) for sec, nsec in zip(data["sec_p"], data["nsec_p"])])
    ts_d = np.array([timestamp(sec, nsec) for sec, nsec in zip(data["sec_d"], data["nsec_d"])])
    ts_diff = np.array([(t_d - t_p).to_sec() * 1e3 for t_p, t_d in zip(ts_p, ts_d)]) # in ms

    pos_p = np.column_stack((data["posx_p"], data["posy_p"], data["posz_p"]))
    pos_d = np.column_stack((data["posx_d"], data["posy_d"], data["posz_d"]))
    distance = np.linalg.norm(pos_p - pos_d, axis=1) / 1000.0
    rho_p = np.linalg.norm(pos_p[:, :2] / 1000.0, axis=1) ** 2
    z_p = pos_p[:, 2] / 1000.0
    rho_d = np.linalg.norm(pos_d[:, :2] / 1000.0, axis=1) ** 2
    z_d = pos_d[:, 2] / 1000.0

    print(f"Loaded {len(data['posx_p'])} events in {filepath}")

    mask_reprod25d = data["run_id"] >= 11049
    data["e_p"][mask_reprod25d] *= 1.00950656406
    data["e_d"][mask_reprod25d] *= 1.00950656406

    e_p_plotter = PromptEnergyPlotter(binmode="nmo")
    e_p_plotter.add(data["e_p"], linecolor="#648fff", fillcolor="#eff3ff")
    e_p_plotter.plot(**meta)
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_e_p_nmo.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_e_p_nmo.png')}")

    e_p_plotter_normal = PromptEnergyPlotter(binmode="normal")
    e_p_plotter_normal.add(data["e_p"], linecolor="#648fff", fillcolor="#eff3ff")
    e_p_plotter_normal.plot(**meta)
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_e_p_normal.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_e_p_normal.png')}")

    e_d_plotter = DelayedEnergyPlotter()
    e_d_plotter.add(data["e_d"], linecolor="#648fff", fillcolor="#eff3ff")
    e_d_plotter.plot()
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_e_d.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_e_d.png')}")

    ts_diff_plotter = PromptDelayedTimePlotter()
    ts_diff_plotter.add(ts_diff, linecolor="#000000", fillcolor="#e5e5e5")
    ts_diff_plotter.plot()
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_dt_p_d.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_dt_p_d.png')}")

    distance_plotter = PromptDelayedDistancePlotter()
    distance_plotter.add(distance, linecolor="#000000", fillcolor="#e5e5e5")
    distance_plotter.plot()
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_dr_p_d.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_dr_p_d.png')}")

    spatial_plotter = SpatialDistributionPlotter()
    spatial_plotter.plot(rho_p, z_p, r"$\rho_{p}$ (m)", r"$z_{p}$ (m)")
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_rho_z_p.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_rho_z_p.png')}")
    spatial_plotter.plot(rho_d, z_d, r"$\rho_{d}$ (m)", r"$z_{d}$ (m)")
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_rho_z_d.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_rho_z_d.png')}")

    xbins = np.linspace(0.0, 1.2, 51)
    ybins = np.linspace(0.0, 3.0, 51)
    muon_veto_plotter = MuonVetoDistributionPlotter()
    muon_veto_plotter.plot(data["dt_mu2p"], data["dlat_mu2p"], xbins, ybins, r"$\Delta t_{\mu-p}$ (s)", r"$d_{\mu-p}$ (m)", is_signal_region=True)
    plt.savefig(f"pdf/{os.path.basename(filepath).replace('.root', '_dt_dlat_p.pdf')}")
    plt.savefig(f"png/{os.path.basename(filepath).replace('.root', '_dt_dlat_p.png')}")

    plt.show()

def analyze_run_info(filepath: str):
    file = uproot.open(filepath)

    tree_daq = file["DAQ"]
    tree_veto = file["Veto"]

    branches_daq = ["run_id", "sec", "nsec"]
    branches_veto = ["run_id", "sec", "nsec", "veto_type", "veto_sec", "veto_nsec"]

    data_daq = tree_daq.arrays(branches_daq, library="np")
    data_veto = tree_veto.arrays(branches_veto, library="np")
    df = tree_veto.arrays(branches_veto, library="pd")

    run_ids = data_daq["run_id"]
    unique_run_ids = np.unique(run_ids)

    min_run = run_ids.min()
    max_run = run_ids.max()
    all_runs = np.arange(min_run, max_run + 1)
    daq_hours = np.zeros_like(all_runs, dtype=float)
    veto_hours = {
        1: np.zeros_like(all_runs, dtype=float), # Beginning of a job
        2: np.zeros_like(all_runs, dtype=float), # Missing headers
        3: np.zeros_like(all_runs, dtype=float), # Big gaps
        4: np.zeros_like(all_runs, dtype=float) # Muons
    }

    # Compute total DAQ time for each RUN

    daq_per_run = {}
    for run_id, sec, nsec in zip(data_daq["run_id"], data_daq["sec"], data_daq["nsec"]):
        ts = timestamp(sec, nsec)
        if run_id not in daq_per_run:
            daq_per_run[run_id] = ts
        else:
            daq_per_run[run_id] += ts

    for run_id in run_ids:
        if run_id in daq_per_run:
            daq_hours[run_id - min_run] = daq_per_run[run_id].to_sec() / 3600.0
        else:
            daq_hours[run_id - min_run] = 0.0

    print(f"Total DAQ time: {np.sum(daq_hours)} hours")

    # Compute total veto time for each veto type

    global_veto_by_type = {}
    global_veto_sums = df.groupby("veto_type")[["veto_sec", "veto_nsec"]].sum()
    for v_type, row in global_veto_sums.iterrows():
        total_ts = timestamp(row["veto_sec"], row["veto_nsec"])
        total_houes = total_ts.to_sec() / 3600.0
        global_veto_by_type[v_type] = total_houes
        print(f"Total veto time for type {v_type}: {total_houes} hours")
    print(f"Total veto time: {np.sum(list(global_veto_by_type.values()))} hours")

    # Compute total veto time for each run and veto type

    aggregated = df.groupby(["run_id", "veto_type"])[["veto_sec", "veto_nsec"]].sum()
    results = {}
    for (run_id, v_type), row in aggregated.iterrows():
        total_duration = timestamp(row["veto_sec"], row["veto_nsec"])
        if run_id not in results:
            results[run_id] = {}
        results[run_id][v_type] = total_duration

    for run, vetos in results.items():
        for v_type, duration in vetos.items():
            veto_hours[v_type][run - min_run] = duration.to_sec() / 3600.0

    fig, ax = plt.subplots(figsize=(16, 6))

    ax.bar(all_runs, daq_hours, width=1.0, align="center", color="#648fff")
    ax.bar(all_runs, veto_hours[3], width=1.0, align="center", color="#ff6464")

    ax.add_patch(
        Rectangle(
            (9789, 0), 11039 - 9789, ax.get_ylim()[1],
            facecolor="#eadcff", edgecolor="none",
            alpha=0.35, zorder=0
        )
    )

    ax.text(
        (9789 + 11039) / 2,
        ax.get_ylim()[1] * 0.92,
        "ReProd25C",
        color="#5b2ca3",
        ha="center",
        va="top",
        fontsize=13,
        fontweight="bold"
    )

    ax.add_patch(
        Rectangle(
            (11049, 0), 12135 - 11049, ax.get_ylim()[1], 
            facecolor="#dbe9ff", edgecolor="none",
            alpha=0.35,
            zorder=0
        )
    )

    ax.text(
        (11049 + 12135) / 2,
        ax.get_ylim()[1] * 0.92,
        "ReProd25D",
        color="#1f4fa3",
        ha="center",
        va="top",
        fontsize=13,
        fontweight="bold"
    )

    ax.set_xlabel("Run ID")
    ax.set_ylabel("Livetime (hours)")

    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(direction="in", which="both", top=True, right=True)

    fig.tight_layout()
    plt.show()
    
    # beg_per_run = {}
    # end_per_run = {}
    # for run_id, sec, nsec in zip(data_veto["run_id"], data_veto["sec"], data_veto["nsec"]):
    #     ts = timestamp(sec, nsec)
    #     if run_id not in beg_per_run:
    #         beg_per_run[run_id] = ts
    #         end_per_run[run_id] = ts
    #     else:
    #         beg_per_run[run_id] = min(beg_per_run[run_id], ts)
    #         end_per_run[run_id] = max(end_per_run[run_id], ts)

if __name__ == "__main__":
    args = parse_args()
    set_latex_style()
    if args.ibd_analysis:
        for filepath in args.ibd_analysis:
            ibd_analysis_plot(filepath) # , reprod="ReProd25D", min_run=11049, max_run=12135)
    if args.cosmo_shape_analysis:
        for filepath in args.cosmo_shape_analysis:
             cosmo_shape_analysis_plot(filepath)
    if args.run_info:
        analyze_run_info(args.run_info)