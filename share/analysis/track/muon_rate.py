import argparse
from datetime import datetime

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
    parser.add_argument("--input", type=str, nargs="+", help="Filepath")
    return parser.parse_args()

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

def calculate_muon_rate(filepath : str):
    file = uproot.open(filepath)
    tree = file["muons"]
    branches = [
        "run_id", "sec", "nsec", # "totq_cd", "totq_wp", 
        "method", "det", "quality",
        "iposx", "iposy", "iposz", 
        "fposx", "fposy", "fposz"
    ]
    data = tree.arrays(branches, library="np")

    mask_cd = np.array([np.any((arr & 1) == 1) for arr in data["det"]])
    mask_wp = np.array([np.any((arr & 2) == 2) for arr in data["det"]])

    mask_cd_wp = np.logical_and(mask_cd, mask_wp)
    mask_cd_only = np.logical_and(mask_cd, np.logical_not(mask_wp))
    mask_wp_only = np.logical_and(np.logical_not(mask_cd), mask_wp)

    data_cd_wp = {key: val[mask_cd_wp] for key, val in data.items()}
    data_cd_only = {key: val[mask_cd_only] for key, val in data.items()}
    data_wp_only = {key: val[mask_wp_only] for key, val in data.items()}

    ts_cd_wp = np.array([timestamp(sec, nsec) for sec, nsec in zip(data_cd_wp["sec"], data_cd_wp["nsec"])])
    ts_cd_only = np.array([timestamp(sec, nsec) for sec, nsec in zip(data_cd_only["sec"], data_cd_only["nsec"])])
    ts_wp_only = np.array([timestamp(sec, nsec) for sec, nsec in zip(data_wp_only["sec"], data_wp_only["nsec"])])

    ts_diff_cd_wp = np.diff(ts_cd_wp)
    ts_diff_cd_only = np.diff(ts_cd_only)
    ts_diff_wp_only = np.diff(ts_wp_only)

    ts_diff_cd_wp = np.array([ts.to_sec() for ts in ts_diff_cd_wp])
    ts_diff_cd_only = np.array([ts.to_sec() for ts in ts_diff_cd_only])
    ts_diff_wp_only = np.array([ts.to_sec() for ts in ts_diff_wp_only])

    fig, ax = plt.subplots(figsize=(7, 6))

    bins = np.linspace(0.0, 2.0, 51)
    centers = 0.5 * (bins[1:] + bins[:-1])
    widths = bins[1:] - bins[:-1]
    hist, _ = np.histogram(ts_diff_cd_wp, bins=bins)
    err = np.sqrt(hist)

    mask = hist > 0
    x_fit = centers[mask]
    y_fit = hist[mask]
    yerr_fit = err[mask]

    A0 = np.max(y_fit)
    lam0 = x_fit[np.argmax(y_fit)]
    p0 = [A0, lam0]

    expdecay = lambda x, A, lam: A * np.exp(-lam * x)
    popt, pcov = curve_fit(expdecay, x_fit, y_fit, p0=p0, sigma=yerr_fit, absolute_sigma=True)
    A, lam = popt
    A_err, lam_err = np.sqrt(np.diag(pcov))

    y_model = expdecay(x_fit, A, lam)
    chisq = np.sum(((y_model - y_fit) / yerr_fit)**2)
    ndf = len(y_fit) - len(popt)
    prob = chi2.sf(chisq, ndf)

    linecolor="#000000"
    fillcolor="#e5e5e5"
    
    ax.fill_between(bins, np.r_[hist, hist[-1]], step="post", color=fillcolor, zorder=1)
    ax.errorbar(centers, hist, yerr=err, xerr=widths/2, fmt="o", color=linecolor, markersize=4.5, zorder=3)

    x_smooth = np.linspace(bins[0], bins[-1], 500)
    y_smooth = expdecay(x_smooth, *popt)

    ax.plot(x_smooth, y_smooth, linestyle="--", linewidth=1.6, color=linecolor, zorder=4)
    text = (
        r"$\chi^2/\mathrm{ndf} = %.1f / %d$" "\n"
        r"$p = %.3f$" "\n\n"
        r"$A = %.2f \pm %.2f$" "\n"
        r"$\lambda = %.2f \pm %.2f~\mathrm{Hz}$"
    ) % (chisq, ndf, prob, A, A_err, lam, lam_err)
    ax.text(0.6, 0.9, text, transform=ax.transAxes, fontsize=15, verticalalignment="top", horizontalalignment="left")

    ax.set_xlabel(r"$\Delta t_{mu}$ (s)")
    ax.set_ylabel(r"Entries")
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(direction="in", which="both", top=True, right=True)
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(bottom=0.0, top=None)

    fig.tight_layout()
    fig.show()

if __name__ == "__main__":
    args = parse_args()
    set_latex_style()
    for filepath in args.input:
        calculate_muon_rate(filepath)