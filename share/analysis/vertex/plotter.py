import argparse
from datetime import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2
import uproot

def set_latex_style():
    mpl.rcParams.update({
        "text.usetex": True, 
        "font.family": "serif", 
        "font.serif": ["Computer Modern Serif"], 
        "mathtext.fontset": "cm", 

        "font.size": 14, 
        "axes.labelsize": 16, 
        "axes.titlesize": 16, 
        "xtick.labelsize": 13, 
        "ytick.labelsize": 13, 
        "legend.fontsize": 13, 

        "axes.linewidth": 1.25, 
        "xtick.direction": "in", 
        "ytick.direction": "in", 
        "xtick.major.size": 8,
        "ytick.major.size": 8,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.minor.width": 0.75,
        "ytick.minor.width": 0.75,
        "xtick.top": True,
        "ytick.right": True,

        "legend.frameon": False,

        "figure.figsize": (8, 6),
        "figure.dpi": 120,

        "savefig.bbox": "tight",
        "savefig.dpi": 300,

    })

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ibd-analysis", type=str, nargs="+", help="IBD analysis filepath")
    parser.add_argument("--cosmo-shape-analysis", type=str, nargs="+", help="Cosmo shape analysis filepath")
    parser.add_argument("--run-info", type=str, help="Run info filepath")
    return parser.parse_args()

def nmo_analysis_bins():
    edges = np.array([0.7, 1.0, 6.6, 7.4, 7.7, 8.1, 8.6, 9.4, 12.0])
    nbins = np.array([1, 56, 4, 1, 1, 1, 1, 1])
    bins = []

    for k in range(len(nbins)):
        start = edges[k]
        end = edges[k + 1]
        nbin = nbins[k]
        segment = np.linspace(start, end, nbin + 1)
        if k > 0:
            segment = segment[1:]
        bins.append(segment)

    return np.concatenate(bins)  

class timestamp:
    def __init__(self, sec : int = 0, nsec : int = 0):
        self.sec = sec
        self.nsec = nsec
        self._normalize()

    def __add__(self, other) -> "timestamp":
        return timestamp(self.sec + other.sec, self.nsec + other.nsec)

    def __sub__(self, other) -> "timestamp":
        return timestamp(self.sec - other.sec, self.nsec - other.nsec)
    
    def __lt__(self, other) -> bool:
        return self.sec < other.sec or (self.sec == other.sec and self.nsec < other.nsec)
    
    def __le__(self, other) -> bool:
        return self.sec < other.sec or (self.sec == other.sec and self.nsec <= other.nsec)
    
    def __gt__(self, other) -> bool:
        return self.sec > other.sec or (self.sec == other.sec and self.nsec > other.nsec)
    
    def __ge__(self, other) -> bool:
        return self.sec > other.sec or (self.sec == other.sec and self.nsec >= other.nsec)

    def __eq__(self, other) -> bool:
        return self.sec == other.sec and self.nsec == other.nsec

    def __ne__(self, other) -> bool:
        return self.sec != other.sec or self.nsec != other.nsec
    
    def __str__(self) -> str:
        dt = datetime.fromtimestamp(self.sec)
        date_str = dt.strftime("%Y-%m-%d %H:%M:%S")
        return f"{date_str}.{self.nsec:09d}"
    
    def to_sec(self) -> float:
        return self.sec + self.nsec / 1e9
    
    def to_nsec(self) -> float:
        return self.sec * 1e9 + self.nsec

    def _normalize(self):
        while self.nsec >= 1e9:
            self.nsec -= 1e9
            self.sec += 1
        while self.nsec < 0:
            self.nsec += 1e9
            self.sec -= 1

class BasePlotter:
    def __init__(self, xlabel="", ylabel="", xlim=None, ylim=(0, None)):
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xlim = xlim
        self.ylim = ylim

    def apply_style(self, ax):
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.tick_params(direction="in", which="both", top=True, right=True)
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)
        if self.xlim: ax.set_xlim(self.xlim)
        if self.ylim: ax.set_ylim(bottom=self.ylim[0], top=self.ylim[1])

class Histogram1DPlotter(BasePlotter):
    def __init__(self, bins, **kwargs):
        super().__init__(**kwargs)
        self.bins = bins
        self.centers = 0.5 * (self.bins[1:] + self.bins[:-1])
        self.widths = self.bins[1:] - self.bins[:-1]
        self.datasets = []

    def add(self, data, linecolor, fillcolor=None, label=None):
        hist, _ = np.histogram(data, bins=self.bins)
        self.datasets.append({
            "hist": hist,
            "err": np.sqrt(hist),
            "linecolor": linecolor,
            "fillcolor": fillcolor,
            "label": label
        })

    def plot(self):
        fig, ax = plt.subplots(figsize=(7, 6))
        
        for d in self.datasets:
            self._draw_dataset(ax, d)
        self._maybe_plot_diff(ax)

        self.apply_style(ax)
        self._maybe_fit(ax)
        
        if any(d["label"] for d in self.datasets):
            ax.legend(loc="upper right")
        
        fig.tight_layout()
        fig.show()

    def _draw_dataset(self, ax, d):
        linestyle = ":" if len(self.datasets) == 2 else "-"
        if d["fillcolor"]:
            ax.fill_between(
                self.bins, np.r_[d["hist"], d["hist"][-1]], 
                step="post", color=d["fillcolor"], zorder=1
            )
        else:
            ax.step(
                self.bins, np.r_[d["hist"], d["hist"][-1]], 
                where="post", color=d["linecolor"], linestyle=linestyle, linewidth=1.2, zorder=2
            )
        
        ax.errorbar(
            self.centers, d["hist"], yerr=d["err"], xerr=self.widths/2, 
            label=d["label"], fmt="o", color=d["linecolor"], markersize=4.5, zorder=3
        )

    def _maybe_fit(self, ax):
        pass

    def _maybe_plot_diff(self, ax):
        if len(self.datasets) != 2: return
        hist1 = self.datasets[0]["hist"]
        err1 = self.datasets[0]["err"]
        hist2 = self.datasets[1]["hist"]
        err2 = self.datasets[1]["err"]

        diff = hist1 - hist2
        err = np.sqrt(err1**2 + err2**2)

        ax.step(
            self.bins, np.r_[diff, diff[-1]], 
            where="post", color="#000000", linewidth=1.2, zorder=2
        )
        ax.errorbar(
            self.centers, diff, yerr=err, xerr=self.widths/2, 
            label="Difference", fmt="o", color="#000000", markersize=4.5, zorder=3
        )

class PromptEnergyPlotter(Histogram1DPlotter):
    def __init__(self, binmode="nmo"):
        bins = nmo_analysis_bins() if binmode == "nmo" else np.linspace(0, 12, 51)
        super().__init__(
            bins=bins, 
            xlabel=r"$E_{p}$ (MeV)", ylabel="Entries", xlim=(0, 12.5)
        )

class DelayedEnergyPlotter(Histogram1DPlotter):
    def __init__(self):
        super().__init__(
            bins=np.linspace(2.0, 2.5, 51), 
            xlabel=r"$E_{d}$ (MeV)", ylabel="Entries", xlim=(1.98, 2.52)
        )

    def _maybe_fit(self, ax):
        if len(self.datasets) != 1: return
        hist = self.datasets[0]["hist"]
        err = self.datasets[0]["err"]
        linecolor = self.datasets[0]["linecolor"]


        mask = hist > 0
        x_fit = self.centers[mask]
        y_fit = hist[mask]
        yerr_fit = err[mask]

        A0 = np.max(y_fit)
        mu0 = x_fit[np.argmax(y_fit)]
        sigma0 = np.std(np.repeat(x_fit, y_fit.astype(int)))
        p0 = [A0, mu0, sigma0]

        gaussian = lambda x, A, mu, sigma: A * np.exp(-(x - mu)**2 / (2 * sigma**2))
        popt, pcov = curve_fit(gaussian, x_fit, y_fit, p0=p0, sigma=yerr_fit, absolute_sigma=True)
        A, mu, sigma = popt
        A_err, mu_err, sigma_err = np.sqrt(np.diag(pcov))

        y_model = gaussian(x_fit, A, mu, sigma)
        chisq = np.sum(((y_model - y_fit) / yerr_fit)**2)
        ndf = len(y_fit) - len(popt)
        prob = chi2.sf(chisq, ndf)

        x_smooth = np.linspace(self.bins[0], self.bins[-1], 500)
        y_smooth = gaussian(x_smooth, *popt)

        ax.plot(
            x_smooth,
            y_smooth,
            linestyle="--",
            linewidth=1.6,
            color=linecolor,
            zorder=4
        )
        text = (
            r"$\chi^2/\mathrm{ndf} = %.1f / %d$" "\n"
            r"$p = %.3f$" "\n\n"
            r"$A = %.2f \pm %.2f$" "\n"
            r"$\mu = %.2f \pm %.2f~\mathrm{MeV}$" "\n"
            r"$\sigma = %.2f \pm %.2f~\mathrm{MeV}$"
        ) % (chisq, ndf, prob, A, A_err, mu, mu_err, sigma, sigma_err)

        ax.text(
            0.6, 0.9,
            text,
            transform=ax.transAxes,
            fontsize=15,
            verticalalignment="top",
            horizontalalignment="left"
        )

class PromptDelayedTimePlotter(Histogram1DPlotter):
    def __init__(self):
        super().__init__(
            bins=np.linspace(0.0, 1.0, 51), 
            xlabel=r"$\Delta t_{p-d}$ (ms)", ylabel="Entries", xlim=(0, 1.02)
        )

    def _maybe_fit(self, ax):
        if len(self.datasets) != 1: return
        hist = self.datasets[0]["hist"]
        err = self.datasets[0]["err"]
        linecolor = self.datasets[0]["linecolor"]

        mask = hist > 0
        x_fit = self.centers[mask]
        y_fit = hist[mask]
        yerr_fit = err[mask]

        mask = x_fit > 0.025
        x_fit = x_fit[mask]
        y_fit = y_fit[mask]
        yerr_fit = yerr_fit[mask]

        A0 = np.max(y_fit)
        tau0 = np.std(np.repeat(x_fit, y_fit.astype(int)))
        p0 = [A0, tau0]

        exponential = lambda x, A, tau: A * np.exp(-x / tau)
        popt, pcov = curve_fit(exponential, x_fit, y_fit, p0=p0, sigma=yerr_fit, absolute_sigma=True)
        A, tau = popt
        A_err, tau_err = np.sqrt(np.diag(pcov))

        y_model = exponential(x_fit, A, tau)
        chisq = np.sum(((y_model - y_fit) / yerr_fit)**2)
        ndf = len(y_fit) - len(popt)
        prob = chi2.sf(chisq, ndf)

        x_smooth = np.linspace(0.025, self.bins[-1], 500)
        y_smooth = exponential(x_smooth, *popt)
    
        ax.plot(
            x_smooth,
            y_smooth,
            linestyle="--",
            linewidth=1.6,
            color="#000000",
            zorder=4
        )
        text = (
            r"$\chi^2/\mathrm{ndf} = %.1f / %d$" "\n"
            r"$p = %.3f$" "\n\n"
            r"$A = %.2f \pm %.2f$" "\n"
            r"$\tau = %.2f \pm %.2f~\mu\mathrm{s}$"
        ) % (chisq, ndf, prob, A, A_err, tau * 1e3, tau_err * 1e3)

        ax.text(
            0.6, 0.9,
            text,
            transform=ax.transAxes,
            fontsize=15,
            verticalalignment="top",
            horizontalalignment="left"
        )

class PromptDelayedDistancePlotter(Histogram1DPlotter):
    def __init__(self):
        super().__init__(
            bins=np.linspace(0.0, 1.5, 51), 
            xlabel=r"$\Delta r_{p-d}$ (m)", ylabel="Entries", xlim=(0, 1.55)
        )

class SpatialDistributionLikePlotter:
    def __init__(self):
        self.xbins = np.linspace(0.0 * 0.0, 17.7 * 17.7, 51)
        self.ybins = np.linspace(-17.7, 17.7, 51)
        self.acrylic_z_curve = np.linspace(-17.7, 17.7, 500)
        self.acrylic_rho2_curve = 17.7**2 - self.acrylic_z_curve**2
        self.fv_z_curve = np.linspace(-16.5, 16.5, 500)
        self.fv_rho2_curve = 16.5**2 - self.fv_z_curve**2

    def plot(self, rho : np.ndarray, z : np.ndarray, xlabel = None, ylabel = None):
        fig, ax = plt.subplots(figsize=(7, 6))
        
        h = ax.hist2d(
            rho,
            z,
            bins=(self.xbins, self.ybins),
            cmin=1.0,
            cmap="jet"
        )
        ax.plot(
            self.acrylic_rho2_curve,
            self.acrylic_z_curve,
            linestyle="--",
            linewidth=1.6,
            color="black"
        )
        ax.plot(
            self.fv_rho2_curve,
            self.fv_z_curve,
            linestyle="--",
            linewidth=1.2,
            color="red"
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        ax.minorticks_on()
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.tick_params(direction="in", which="both", top=True, right=True)
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)

        cbar = fig.colorbar(h[3], ax=ax)
        cbar.set_label(r"Entries")
        cbar.ax.minorticks_on()
        cbar.ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        cbar.ax.tick_params(
            direction="in",
            which="both",
            width=1.2,
        )

        ax.set_xlim(0.0, 320.0)
        ax.set_ylim(-17.8, 17.8)

        fig.tight_layout()
        fig.show()

def cosmo_shape_analysis_plot(filepath: str, **meta):
    file = uproot.open(filepath)
    tree_bkg = file["background_events"]
    tree_sig = file["signal_events"]

    branches = [
        "posx_p", "posy_p", "posz_p", "sec_p", "nsec_p", "e_p", "dlat_mu2p", "dt_mu2p",
        "posx_d", "posy_d", "posz_d", "sec_d", "nsec_d", "e_d", "dlat_mu2d", "dt_mu2d"
    ]

    data_bkg = tree_bkg.arrays(branches, library="np")
    data_sig = tree_sig.arrays(branches, library="np")

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

    e_p_plotter = PromptEnergyPlotter(binmode="normal")
    e_p_plotter.add(data_sig["e_p"], linecolor="#648fff", label="Cosmogenic enriched region")
    e_p_plotter.add(data_bkg["e_p"], linecolor="#ff6464", label="Cosmogenic depleted region")
    e_p_plotter.plot(**meta)

    e_d_plotter = DelayedEnergyPlotter()
    e_d_plotter.add(data_sig["e_d"], linecolor="#648fff", label="Cosmogenic enriched region")
    e_d_plotter.add(data_bkg["e_d"], linecolor="#ff6464", label="Cosmogenic depleted region")
    e_d_plotter.plot()

    ts_diff_plotter = PromptDelayedTimePlotter()
    ts_diff_plotter.add(ts_diff_sig, linecolor="#648fff", label="Cosmogenic enriched region")
    ts_diff_plotter.add(ts_diff_bkg, linecolor="#ff6464", label="Cosmogenic depleted region")
    ts_diff_plotter.plot()

    distance_plotter = PromptDelayedDistancePlotter()
    distance_plotter.add(distance_sig, linecolor="#648fff", label="Cosmogenic enriched region")
    distance_plotter.add(distance_bkg, linecolor="#ff6464", label="Cosmogenic depleted region")
    distance_plotter.plot()

    spatial_plotter = SpatialDistributionLikePlotter()
    spatial_plotter.plot(rho_p_bkg, z_p_bkg, r"$\rho_{p}$ (m)", r"$z_{p}$ (m)")
    spatial_plotter.plot(rho_d_bkg, z_d_bkg, r"$\rho_{d}$ (m)", r"$z_{d}$ (m)")
    spatial_plotter.plot(rho_p_sig, z_p_sig, r"$\rho_{p}$ (m)", r"$z_{p}$ (m)")
    spatial_plotter.plot(rho_d_sig, z_d_sig, r"$\rho_{d}$ (m)", r"$z_{d}$ (m)")

    plt.show()

def ibd_analysis_plot(filepath: str, **meta):
    file = uproot.open(filepath)
    tree = file["events"]

    branches = [
        "posx_p", "posy_p", "posz_p", "sec_p", "nsec_p", "e_p",
        "posx_d", "posy_d", "posz_d", "sec_d", "nsec_d", "e_d"
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

    e_p_plotter = PromptEnergyPlotter(binmode="nmo")
    e_p_plotter.add(data["e_p"], linecolor="#648fff", fillcolor="#eff3ff")
    e_p_plotter.plot(**meta)

    e_d_plotter = DelayedEnergyPlotter()
    e_d_plotter.add(data["e_d"], linecolor="#648fff", fillcolor="#eff3ff")
    e_d_plotter.plot()

    ts_diff_plotter = PromptDelayedTimePlotter()
    ts_diff_plotter.add(ts_diff, linecolor="#000000", fillcolor="#e5e5e5")
    ts_diff_plotter.plot()

    distance_plotter = PromptDelayedDistancePlotter()
    distance_plotter.add(distance, linecolor="#000000", fillcolor="#e5e5e5")
    distance_plotter.plot()

    spatial_plotter = SpatialDistributionLikePlotter()
    spatial_plotter.plot(rho_p, z_p, r"$\rho_{p}$ (m)", r"$z_{p}$ (m)")
    spatial_plotter.plot(rho_d, z_d, r"$\rho_{d}$ (m)", r"$z_{d}$ (m)")

    plt.show()

def ibd_rate_plot(filepath: str, run_info_filepath: str):
    file = uproot.open(filepath)
    tree = file["events"]

    branches = [
        "posx_p", "posy_p", "posz_p", "sec_p", "nsec_p", "e_p",
        "posx_d", "posy_d", "posz_d", "sec_d", "nsec_d", "e_d"
    ]
    data = tree.arrays(branches, library="np")

    file = uproot.open(run_info_filepath)
    tree = file["run_info"]
    branches = ["run_id", "earliest_sec", "earliest_nsec", "latest_sec", "latest_nsec"]
    run_info = tree.arrays(branches, library="np")
    earliest_ts = np.array([timestamp(sec, nsec) for sec, nsec in zip(run_info["earliest_sec"], run_info["earliest_nsec"])])
    latest_ts = np.array([timestamp(sec, nsec) for sec, nsec in zip(run_info["latest_sec"], run_info["latest_nsec"])])
    duration = latest_ts - earliest_ts

    run_id = []
    for i in range(len(data["posx_p"])):
        for j in range(len(earliest_ts)):
            if earliest_ts[j] <= timestamp(data["sec_p"][i], data["nsec_p"][i]) <= latest_ts[j]:
                run_id.append(run_info["run_id"][j])
                break
    
    data["run_id"] = np.array(run_id)

    per_run_nevents = []
    for run in run_info["run_id"]:
        mask = data["run_id"] == run
        ts_p = np.array([timestamp(sec, nsec) for sec, nsec in zip(data["sec_p"][mask], data["nsec_p"][mask])])
        per_run_nevents.append(len(ts_p))

    per_run_ibd_rate = []
    per_run_ibd_rate_err = []
    for run, nevents, dur in zip(run_info["run_id"], per_run_nevents, duration):
        ibd_rate = nevents / (dur.to_sec() / 3600.0 / 24.0)
        per_run_ibd_rate.append(ibd_rate)
        per_run_ibd_rate_err.append(np.sqrt(nevents) / (dur.to_sec() / 3600.0 / 24.0))

    Ntot = np.sum(per_run_nevents)
    duration_sec = np.array([d.to_sec() for d in duration])
    durtot = np.sum(duration_sec) / 3600.0 / 24.0
    print(f"Total ({Ntot} events): rate = {Ntot / durtot:.3f} +- {np.sqrt(Ntot) / durtot:.3f} events/day")

    # for run, nevents, rate, err in zip(run_info["run_id"], per_run_nevents, per_run_ibd_rate, per_run_ibd_rate_err):
    #     print(f"RUN {run} ({nevents} events): rate = {rate:.3f} +- {err:.3f} events/day")

    plt.figure()
    plt.errorbar(run_info["run_id"], per_run_ibd_rate, yerr=per_run_ibd_rate_err, fmt="o")
    plt.show()

if __name__ == "__main__":
    args = parse_args()
    set_latex_style()
    # if args.ibd_analysis:
    #     for filepath in args.ibd_analysis:
    #         ibd_analysis_plot(filepath) # , reprod="ReProd25D", min_run=11049, max_run=12135)
    # if args.cosmo_shape_analysis:
    #     for filepath in args.cosmo_shape_analysis:
    #         cosmo_shape_analysis_plot(filepath)
    if args.run_info:
        ibd_rate_plot(args.ibd_analysis[0], args.run_info)