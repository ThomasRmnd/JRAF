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
        "text.usetex": False, # TODO: need to be changed for good plots
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
    parser.add_argument("--run", type=int, help="Run number")
    # parser.add_argument("--input", type=str, nargs="+", help="Filepath")
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

def exponential_decay(x, A, lam):
    return A * np.exp(-lam * x)

def fit_exponential_decay(x, y, yerr):
    A0 = np.max(y)
    lam0 = np.std(np.repeat(x, y.astype(int)))
    p0 = [A0, lam0]

    popt, pcov = curve_fit(exponential_decay, x, y, p0=p0, sigma=yerr, absolute_sigma=True)
    A, lam = popt
    A_err, lam_err = np.sqrt(np.diag(pcov))

    y_model = exponential_decay(x, A, lam)
    chisq = np.sum(((y_model - y) / yerr)**2)
    ndf = len(y) - len(popt)
    prob = chi2.sf(chisq, ndf)

    return chisq, ndf, prob, A, lam, A_err, lam_err

# def calculate_muon_rate(filepath : str, plot=False):
def calculate_muon_rate(run : int, plot=False):
    # file = uproot.open(filepath)
    file = uproot.open(f"/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/summary/RUN.{run}.output.reprod25c.cca.root")
    tree = file["muons"]
    branches = ["run_id", "sec", "nsec", "totq_cd", "totq_wp"]
    data = tree.arrays(branches, library="np")

    # mask_cd = np.array([np.any((arr & 1) == 1) for arr in data["det"]])
    # mask_wp = np.array([np.any((arr & 2) == 2) for arr in data["det"]])

    mask_cd_wp = np.logical_and(data["totq_cd"] > 0, data["totq_wp"] > 0)
    mask_cd_only = np.logical_and(data["totq_cd"] > 0, data["totq_wp"] == 0)
    mask_wp_only = np.logical_and(data["totq_cd"] == 0, data["totq_wp"] > 0)

    data_cd_wp = {key: val[mask_cd_wp] for key, val in data.items()}
    data_cd_only = {key: val[mask_cd_only] for key, val in data.items()}
    data_wp_only = {key: val[mask_wp_only] for key, val in data.items()}

    ts_cd_wp = np.array([timestamp(sec, nsec) for sec, nsec in zip(data_cd_wp["sec"], data_cd_wp["nsec"])])
    ts_cd_only = np.array([timestamp(sec, nsec) for sec, nsec in zip(data_cd_only["sec"], data_cd_only["nsec"])])
    ts_wp_only = np.array([timestamp(sec, nsec) for sec, nsec in zip(data_wp_only["sec"], data_wp_only["nsec"])])

    ts_min = ts_cd_wp.min()
    ts_max = ts_cd_wp.max()
    if ts_cd_only.size > 0:
        ts_min = min(ts_min, ts_cd_only.min())
        ts_max = max(ts_max, ts_cd_only.max())
    if ts_wp_only.size > 0:
        ts_min = min(ts_min, ts_wp_only.min())
        ts_max = max(ts_max, ts_wp_only.max())

    ts_diff_cd_wp = np.diff(ts_cd_wp)
    ts_diff_cd_only = np.diff(ts_cd_only)
    ts_diff_wp_only = np.diff(ts_wp_only)

    ts_diff_cd_wp = np.array([ts.to_sec() for ts in ts_diff_cd_wp])
    ts_diff_cd_only = np.array([ts.to_sec() for ts in ts_diff_cd_only])
    ts_diff_wp_only = np.array([ts.to_sec() for ts in ts_diff_wp_only])

    bins = np.linspace(0.0, 2.0, 51)
    centers = 0.5 * (bins[1:] + bins[:-1])
    widths = bins[1:] - bins[:-1]
    
    hist_cd_wp, _ = np.histogram(ts_diff_cd_wp, bins=bins)
    err_cd_wp = np.sqrt(hist_cd_wp)
    hist_cd_only, _ = np.histogram(ts_diff_cd_only, bins=bins)
    err_cd_only = np.sqrt(hist_cd_only)
    hist_wp_only, _ = np.histogram(ts_diff_wp_only, bins=bins)
    err_wp_only = np.sqrt(hist_wp_only)

    ts_diffs = [ts_diff_cd_wp, ts_diff_cd_only, ts_diff_wp_only]
    hists = [hist_cd_wp, hist_cd_only, hist_wp_only]
    errs = [err_cd_wp, err_cd_only, err_wp_only]
    rates = []
    rates_err = []
    for diff, h, e in zip(ts_diffs, hists, errs):
        mask = h > 0
        x_fit = centers[mask]
        y_fit = h[mask]
        yerr_fit = e[mask]

        if len(x_fit) < 2:
            # number of counts divided by (ts_max - ts_min).to_sec()
            rates.append(len(diff) / (ts_max - ts_min).to_sec())
            rates_err.append(np.sqrt(len(diff)) / (ts_max - ts_min).to_sec())
        else:
            chisq, ndf, prob, A, lam, A_err, lam_err = fit_exponential_decay(x_fit, y_fit, yerr_fit)
            rates.append(lam)
            rates_err.append(lam_err)

    print(f"Run ID: {data['run_id'][0]}")
    print(f"Timestamp: {datetime.fromtimestamp(data['sec'][0]).strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"CD-WP rate: {rates[0]:.2f} +/- {rates_err[0]:.2f} cps")
    print(f"CD only rate: {rates[1]:.2f} +/- {rates_err[1]:.2f} cps")
    print(f"WP only rate: {rates[2]:.2f} +/- {rates_err[2]:.2f} cps")

    # Save the rates in file
    ofile=f"/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/summary/rates/RUN.{run}.rates.root"
    with uproot.recreate(ofile) as f:
        f["rates"] = {
            "run_id": np.full(len(rates), data["run_id"][0], dtype=np.int32),
            "sec": np.full(len(rates), data["sec"][0], dtype=np.int64),
            "rate_cd_wp": rates[0].astype(np.float64),
            "rate_cd_wp_err": rates_err[0].astype(np.float64),
            "rate_cd_only": rates[1].astype(np.float64),
            "rate_cd_only_err": rates_err[1].astype(np.float64),
            "rate_wp_only": rates[2].astype(np.float64),
            "rate_wp_only_err": rates_err[2].astype(np.float64)
        }

    if not plot:
        return data["run_id"][0], data["sec"][0], rates, rates_err

    linecolors = ["#000000", "#648fff", "#ff6464"]
    fillcolors = ["#e5e5e5", "#eff3ff", "#ffefef"]

    fig, ax = plt.subplots(figsize=(7, 6))

    for h, e, lcolor, fcolor in zip(hists, errs, linecolors, fillcolors):
        ax.fill_between(bins, np.r_[h, h[-1]], step="post", color=fcolor, zorder=1)
        ax.errorbar(centers, h, yerr=e, xerr=widths/2, fmt="o", color=lcolor, markersize=4.5, zorder=3)

        mask = h > 0
        x_fit = centers[mask]
        y_fit = h[mask]
        yerr_fit = e[mask]

        if len(x_fit) < 2:
            continue

        chisq, ndf, prob, A, lam, A_err, lam_err = fit_exponential_decay(x_fit, y_fit, yerr_fit)

        x_smooth = np.linspace(bins[0], bins[-1], 500)
        y_smooth = exponential_decay(x_smooth, A, lam)

        ax.plot(x_smooth, y_smooth, linestyle="--", linewidth=1.6, color=lcolor, zorder=4)
        # text = (
        #     r"$\chi^2/\mathrm{ndf} = %.1f / %d$" "\n"
        #     r"$p = %.3f$" "\n\n"
        #     r"$A = %.2f \pm %.2f$" "\n"
        #     r"$\lambda = %.2f \pm %.2f~\mathrm{Hz}$"
        # ) % (chisq, ndf, prob, A, A_err, lam, lam_err)
        # ax.text(0.6, 0.9, text, transform=ax.transAxes, fontsize=15, verticalalignment="top", horizontalalignment="left")

    ax.set_xlabel(r"$\Delta t_{\mu}$ (s)")
    ax.set_ylabel(r"Entries")
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(direction="in", which="both", top=True, right=True)
    ax.set_xlim(0.0, 2.0)
    ax.set_yscale("log")

    fig.tight_layout()
    plt.show()

    return data["run_id"][0], data["sec"][0], rates, rates_err

if __name__ == "__main__":
    args = parse_args()
    set_latex_style()

    calculate_muon_rate(args.run, plot=True)
    
    '''
    run_ids = []
    timestamps = []
    rates_cd_wp = []
    rates_cd_only = []
    rates_wp_only = []
    rates_err_cd_wp = []
    rates_err_cd_only = []
    rates_err_wp_only = []
    
    for filepath in args.input:
        run_id, ts, rates, rates_err = calculate_muon_rate(filepath, plot=False)
        run_ids.append(run_id)
        timestamps.append(ts)
        rates_cd_wp.append(rates[0])
        rates_cd_only.append(rates[1])
        rates_wp_only.append(rates[2])
        rates_err_cd_wp.append(rates_err[0])
        rates_err_cd_only.append(rates_err[1])
        rates_err_wp_only.append(rates_err[2])

    rates_total = np.array(rates_cd_wp) + np.array(rates_cd_only) + np.array(rates_wp_only)
    err_total = np.sqrt(np.array(rates_err_cd_wp)**2 + np.array(rates_err_cd_only)**2 + np.array(rates_err_wp_only)**2)

    config = [
        (rates_total, err_total, "#000000", f"Total: {np.mean(rates_total):.2f} +/- {np.std(rates_total):.2f} Hz"),
        (rates_cd_wp, rates_err_cd_wp, "#e69f00", f"CD+WP: {np.mean(rates_cd_wp):.2f} +/- {np.std(rates_cd_wp):.2f} Hz"),
        (rates_cd_only, rates_err_cd_only, "#009e73", f"CD only: {np.mean(rates_cd_only):.2f} +/- {np.std(rates_cd_only):.2f} Hz"),
        (rates_wp_only, rates_err_wp_only, "#56b4e9", f"WP only: {np.mean(rates_wp_only):.2f} +/- {np.std(rates_wp_only):.2f} Hz"),
    ]

    fig, ax = plt.subplots(figsize=(16, 6))

    for data, err, color, label in config:
        ax.errorbar(run_ids, data, yerr=err, fmt="o", color=color, markersize=5.0, capsize=0, elinewidth=1.0, markeredgecolor="k", markeredgewidth=0.5, label=label, zorder=3)
        mean_val = np.mean(data)
        ax.axhline(mean_val, color=color, linestyle="--", linewidth=2.0, zorder=2)

    ax.set_xlabel("Run Number", fontsize=14)
    ax.set_ylabel("Muon Rate [Hz]", fontsize=14)
    ax.set_ylim(-0.1, 10)

    ax.tick_params(direction='in', which='both', top=True, right=True, labelsize=12)
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    title_str = f"Run range: {min(run_ids)} - {max(run_ids)}" # exposure: {exposure_days:.1f} days"
    ax.set_title(title_str, loc='right', fontsize=14, color='grey', pad=20)

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=4, frameon=False, fontsize=12, handletextpad=0.1)

    fig.tight_layout()
    plt.show()
    '''