import argparse
from datetime import datetime
import os

import matplotlib as mpl
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chi2
import uproot

def set_latex_style():
    mpl.rcParams.update({
        "text.usetex": True, 
        "font.family": "serif", 
        "font.serif": ["Computer Modern Serif"], 
        "mathtext.fontset": "cm", 

        "font.size": 22, 
        "axes.labelsize": 22, 
        "axes.titlesize": 22, 
        "xtick.labelsize": 18, 
        "ytick.labelsize": 18, 
        "legend.fontsize": 18, 

        "axes.linewidth": 1.35, 
        "xtick.direction": "in", 
        "ytick.direction": "in", 
        "xtick.major.size": 10,
        "ytick.major.size": 10,
        "xtick.minor.size": 5,
        "ytick.minor.size": 5,
        "xtick.major.width": 1.25,
        "ytick.major.width": 1.25,
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
    parser.add_argument("--input", type=str, nargs=2, help="IBD analysis filepath")
    parser.add_argument("--label", type=str, nargs=2, help="Labels for the two datasets")
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
        carry = self.nsec // 1e9
        self.sec += carry
        self.nsec %= 1e9

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
        if len(self.datasets) >= 2:
            print(f"Warning: This plotter is designed for exactly 2 datasets. Ignoring extra dataset")
            return

        hist, _ = np.histogram(data, bins=self.bins)
        self.datasets.append({
            "hist": hist,
            "err": np.sqrt(hist),
            "linecolor": linecolor,
            "fillcolor": fillcolor,
            "label": label
        })

    def plot(self):
        if len(self.datasets) != 2:
            raise ValueError("Plotting requires exactly 2 datasets")

        fig, (ax_main, ax_diff) = plt.subplots(
            2, 1, figsize=(7, 10), sharex=True,
            gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05}
        )
        
        for d in self.datasets:
            self._draw_main(ax_main, d)
        
        self._draw_diff(ax_diff)
        self._apply_shared_style(ax_main, ax_diff)
        
        if any(d["label"] for d in self.datasets):
            ax_main.legend(loc="upper right", frameon=False)
        
        fig.show()

    def _draw_main(self, ax, d):
        if d["fillcolor"]:
            ax.fill_between(
                self.bins, np.r_[d["hist"], d["hist"][-1]], 
                step="post", color=d["fillcolor"], alpha=0.3, zorder=1
            )
        else:
            ax.step(
                self.bins, np.r_[d["hist"], d["hist"][-1]], 
                where="post", color=d["linecolor"], linestyle=":", linewidth=1.2, zorder=2
            )
        
        ax.errorbar(
            self.centers, d["hist"], yerr=d["err"], xerr=self.widths/2, 
            label=d["label"], fmt="o", color=d["linecolor"], markersize=4.5, zorder=3
        )

    def _draw_diff(self, ax):
        h1, e1 = self.datasets[0]["hist"], self.datasets[0]["err"]
        h2, e2 = self.datasets[1]["hist"], self.datasets[1]["err"]

        diff = h1 - h2
        err = np.sqrt(e1**2 + e2**2)

        ax.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.7)
        ax.plot(self.centers, diff, linestyle='None', marker="o", color="black", markersize=4, zorder=3)
        ax.errorbar(self.centers, diff, yerr=err, xerr=self.widths/2, fmt="o", color="black", markersize=4.5, zorder=3)

    def _apply_shared_style(self, ax_main, ax_diff):
        for ax in [ax_main, ax_diff]:
            self.apply_style(ax)
        ax_main.set_xlabel(None)

        ax_main.set_ylabel("Entries")
        ax_diff.set_ylabel(r"$\Delta$")

        ax_diff.yaxis.set_major_locator(plt.MaxNLocator(5, prune='both'))

        h1, e1 = self.datasets[0]["hist"], self.datasets[0]["err"]
        h2, e2 = self.datasets[1]["hist"], self.datasets[1]["err"]
        diff = h1 - h2
        max_deviation = np.max(np.abs(diff))
        limit = max_deviation * 1.25
        ax_diff.set_ylim(bottom=-limit, top=limit)
        
        plt.setp(ax_main.get_xticklabels(), visible=False)

class PromptEnergyPlotter(Histogram1DPlotter):
    def __init__(self, binmode="nmo", **kwargs):
        bins = nmo_analysis_bins() if binmode == "nmo" else np.linspace(0, 12, 51)
        super().__init__(
            bins=bins, 
            xlabel=r"$E_{p}$ (MeV)", ylabel="Entries", xlim=(0, 12.5),
            **kwargs
        )

class DelayedEnergyPlotter(Histogram1DPlotter):
    def __init__(self, **kwargs):
        super().__init__(
            bins=np.linspace(2.0, 2.5, 76), 
            xlabel=r"$E_{d}$ (MeV)", ylabel="Entries", xlim=(1.98, 2.52),
            **kwargs
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
            0.55, 0.9,
            text,
            transform=ax.transAxes,
            fontsize=18,
            verticalalignment="top",
            horizontalalignment="left"
        )

class PromptDelayedTimePlotter(Histogram1DPlotter):
    def __init__(self, **kwargs):
        super().__init__(
            bins=np.linspace(0.0, 1.0, 51), 
            xlabel=r"$\Delta t_{p-d}$ (ms)", ylabel="Entries", xlim=(0, 1.02),
            **kwargs
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
            0.55, 0.9,
            text,
            transform=ax.transAxes,
            fontsize=18,
            verticalalignment="top",
            horizontalalignment="left"
        )

class PromptDelayedDistancePlotter(Histogram1DPlotter):
    def __init__(self, **kwargs):
        super().__init__(
            bins=np.linspace(0.0, 1.5, 51), 
            xlabel=r"$\Delta r_{p-d}$ (m)", ylabel="Entries", xlim=(0, 1.55),
            **kwargs
        )

class SpatialDistributionPlotter:
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

def plot_comparator(filepath1 : str, filepath2 : str, label1 : str, label2 : str):
    file1 = uproot.open(filepath1) 
    file2 = uproot.open(filepath2) 
    
    tree1 = file1["events"]
    tree2 = file2["events"]

    if "IBD_all_reprod" in filepath1:
        branches1 = [
            "run_number",
            "vertex_x_p_omilrec", "vertex_y_p_omilrec", "vertex_z_p_omilrec", "time_p_ns", "energy_p_omilrec",
            "vertex_x_d_omilrec", "vertex_y_d_omilrec", "vertex_z_d_omilrec", "time_d_ns", "energy_d_omilrec"
        ]
    else:
        branches1 = [
            "run_id",
            "posx_p", "posy_p", "posz_p", "sec_p", "nsec_p", "e_p",
            "posx_d", "posy_d", "posz_d", "sec_d", "nsec_d", "e_d"
        ]
    
    if "IBD_all_reprod" in filepath2:
        branches2 = [
            "run_number",
            "vertex_x_p_omilrec", "vertex_y_p_omilrec", "vertex_z_p_omilrec", "time_p_ns", "energy_p_omilrec",
            "vertex_x_d_omilrec", "vertex_y_d_omilrec", "vertex_z_d_omilrec", "time_d_ns", "energy_d_omilrec"
        ]
    else:
        branches2 = [
            "run_id",
            "posx_p", "posy_p", "posz_p", "sec_p", "nsec_p", "e_p",
            "posx_d", "posy_d", "posz_d", "sec_d", "nsec_d", "e_d"
        ]

    data1 = tree1.arrays(branches1, library="np")
    data2 = tree2.arrays(branches2, library="np")

    if "IBD_all_reprod" in filepath1:
        data1["run_number"] = np.round(data1["run_number"]).astype(int)
    if "IBD_all_reprod" in filepath2: 
        data2["run_number"] = np.round(data2["run_number"]).astype(int)

    max_run_id = 11039
    failed_jobs = [
        9803, 10155, 10212, 10231, 10241, 10252, 10261, 10270, 10370, 10390, 
        10459, 10470, 10479, 10520, 10529, 10540, 10550, 10563, 10584, 10593, 
        11027, 11237, 11410, 11397, 11788
    ]

    mask1 = (data1["run_id"] if "IBD_all_reprod" not in filepath1 else data1["run_number"]) <= max_run_id 
    mask2 = (data2["run_id"] if "IBD_all_reprod" not in filepath2 else data2["run_number"]) <= max_run_id 
    data1 = {k: v[mask1] for k, v in data1.items()} 
    data2 = {k: v[mask2] for k, v in data2.items()}

    mask1 = ~np.isin(data1["run_id"] if "IBD_all_reprod" not in filepath1 else data1["run_number"], failed_jobs) 
    mask2 = ~np.isin(data2["run_id"] if "IBD_all_reprod" not in filepath2 else data2["run_number"], failed_jobs)
    data1 = {k: v[mask1] for k, v in data1.items()} 
    data2 = {k: v[mask2] for k, v in data2.items()}

    e_p_plotter = PromptEnergyPlotter(binmode="normal")
    e_p_plotter.add(data1["e_p"] if "IBD_all_reprod" not in filepath1 else data1["energy_p_omilrec"], linecolor="#648fff", fillcolor="#eff3ff", label=label1)
    e_p_plotter.add(data2["e_p"] if "IBD_all_reprod" not in filepath2 else data2["energy_p_omilrec"], linecolor="#ff6464", fillcolor="#ffefef", label=label2)
    e_p_plotter.plot()

    e_d_plotter = DelayedEnergyPlotter()
    e_d_plotter.add(data1["e_d"] if "IBD_all_reprod" not in filepath1 else data1["energy_d_omilrec"], linecolor="#648fff", label=label1)
    e_d_plotter.add(data2["e_d"] if "IBD_all_reprod" not in filepath2 else data2["energy_d_omilrec"], linecolor="#ff6464", label=label2)
    e_d_plotter.plot()

    plt.show()

if __name__ == "__main__":
    args = parse_args()
    set_latex_style()
    plot_comparator(args.input[0], args.input[1], args.label[0], args.label[1])
