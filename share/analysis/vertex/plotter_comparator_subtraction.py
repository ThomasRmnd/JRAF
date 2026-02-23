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

    def add(self, data_sig, data_bkg, linecolor, fillcolor=None, label=None):
        if len(self.datasets) >= 2:
            print(f"Warning: This plotter is designed for exactly 2 datasets. Ignoring extra dataset")
            return

        hist_sig, _ = np.histogram(data_sig, bins=self.bins)
        hist_bkg, _ = np.histogram(data_bkg, bins=self.bins)
        hist = hist_sig - hist_bkg
        err = np.sqrt(hist_sig + hist_bkg)
        self.datasets.append({
            "hist": hist,
            "err": err,
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

def plot_comparator(filepath1 : str, filepath2 : str, label1 : str, label2 : str):
    file1 = uproot.open(filepath1) 
    file2 = uproot.open(filepath2) 
    
    tree1_sig = file1["sig"]
    tree1_bkg = file1["bkg"]
    tree2_sig = file2["sig"]
    tree2_bkg = file2["bkg"]

    branches = [
        "run_id",
        "posx_p", "posy_p", "posz_p", "sec_p", "nsec_p", "e_p",
        "posx_d", "posy_d", "posz_d", "sec_d", "nsec_d", "e_d"
    ]

    data1_sig = tree1_sig.arrays(branches, library="np")
    data1_bkg = tree1_bkg.arrays(branches, library="np")
    data2_sig = tree2_sig.arrays(branches, library="np")
    data2_bkg = tree2_bkg.arrays(branches, library="np")

    e_p_plotter = PromptEnergyPlotter(binmode="normal")
    e_p_plotter.add(data1_sig["e_p"], data1_bkg["e_p"], linecolor="#648fff", fillcolor="#eff3ff", label=label1)
    e_p_plotter.add(data2_sig["e_p"], data2_bkg["e_p"], linecolor="#ff6464", fillcolor="#ffefef", label=label2)
    e_p_plotter.plot()

    plt.show()

if __name__ == "__main__":
    args = parse_args()
    set_latex_style()
    plot_comparator(args.input[0], args.input[1], args.label[0], args.label[1])
