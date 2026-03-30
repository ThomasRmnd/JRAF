import argparse
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FuncFormatter
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chi2, norm
import uproot

def set_latex_style():
    mpl.rcParams.update({
        "text.usetex": True, # TODO: need to be changed for good plots
        "font.family": "serif", 
        "font.serif": ["Computer Modern Serif"], 
        "mathtext.fontset": "cm", 

        "font.size": 20, 
        "axes.labelsize": 20, 
        "axes.titlesize": 18, 
        "xtick.labelsize": 18, 
        "ytick.labelsize": 18, 
        "legend.fontsize": 17, 

        "axes.linewidth": 1.25, 
        "xtick.direction": "in", 
        "ytick.direction": "in", 
        "xtick.major.size": 10,
        "ytick.major.size": 10,
        "xtick.minor.size": 5,
        "ytick.minor.size": 5,
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
    parser.add_argument("--run", type=int, default=0, help="Run number")
    parser.add_argument("--input", type=str, default="", help="Filepath")
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

def find_run_file(run: int) -> Path:
    base_dir = Path("/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/summary")
    
    pattern = f"RUN.{run}.output.reprod*.cca.root"
    matches = list(base_dir.glob(pattern))

    if len(matches) == 0:
        raise FileNotFoundError(f"No file found for run {run} with pattern {pattern}")
    
    if len(matches) > 1:
        raise RuntimeError(
            f"Multiple files found for run {run}:\n" +
            "\n".join(str(m) for m in matches)
        )

    return matches[0]

def apply_veto(ts_target : np.ndarray, ts_all : np.ndarray, veto_window : timestamp = timestamp(0, 10000000)) -> np.ndarray:
    mask_keep = np.ones(len(ts_target), dtype=bool)
    veto_window_full = np.full(len(ts_all), veto_window)
    for i, ts in enumerate(ts_target):
        dt = ts_all - ts
        if np.any((veto_window_full < dt) & (dt < veto_window_full) & (dt != timestamp())):
            mask_keep[i] = False
    return ts_target[mask_keep]

def calculate_wp_tagging_efficiency(ts_all : np.ndarray, ts_cd_only : np.ndarray, rate_total_cd : float, nb_cd_wp_high : int, veto_window : timestamp = timestamp(0, 10000000)):
    prob_veto = 1.0 - np.exp(-rate_total_cd * veto_window.to_sec())
    ts_cd_only_corr = apply_veto(ts_cd_only, ts_all, veto_window)
    nb_cd_only_corr = len(ts_cd_only_corr) * (1.0 - prob_veto)
    return 1.0 - (nb_cd_only_corr / nb_cd_wp_high)

def calculate_muon_rate(run : int, plot=False):
    filepath = find_run_file(run)
    print(f"Using file: {filepath}")

    file = uproot.open(filepath)
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

    names = ["CD+WP", "CD only", "WP only"]
    ts_diffs = [ts_diff_cd_wp, ts_diff_cd_only, ts_diff_wp_only]
    hists = [hist_cd_wp, hist_cd_only, hist_wp_only]
    errs = [err_cd_wp, err_cd_only, err_wp_only]
    rates = []
    rates_err = []
    for name, diff, h, e in zip(names, ts_diffs, hists, errs):
        mask = h > 0
        x_fit = centers[mask]
        y_fit = h[mask]
        yerr_fit = e[mask]

        if name == "CD only":
            duration = (ts_max - ts_min).to_sec()
            rates.append(len(diff) / duration)
            rates_err.append(np.sqrt(len(diff)) / duration)
        else:
            chisq, ndf, prob, A, lam, A_err, lam_err = fit_exponential_decay(x_fit, y_fit, yerr_fit)
            rates.append(lam)
            rates_err.append(lam_err)

    ts_all = np.concatenate([ts_cd_wp, ts_cd_only, ts_wp_only])
    mask_high_charge = data["totq_cd"] > 30000
    mask_cd_wp_high = np.logical_and(mask_cd_wp, mask_high_charge)
    wp_tagging_efficiency = calculate_wp_tagging_efficiency(ts_all, ts_cd_only, rates[0] + rates[1], np.sum(mask_cd_wp_high))

    print(f"Run ID: {data['run_id'][0]}")
    print(f"Timestamp: {datetime.fromtimestamp(data['sec'][0]).strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"CD-WP rate: {rates[0]:.2f} +/- {rates_err[0]:.2f} cps")
    print(f"CD only rate: {rates[1]:.2f} +/- {rates_err[1]:.2f} cps")
    print(f"WP only rate: {rates[2]:.2f} +/- {rates_err[2]:.2f} cps")
    print(f"WP tagging efficiency: {wp_tagging_efficiency * 100.0:.2f}")

    # Save the rates in file
    ofile=f"/sps/juno/jdeandre/rtraw_ThomasRaymond/reconstruction/reprod/summary/rates/RUN.{run}.rates.root"
    with uproot.recreate(ofile) as f:
        f["rates"] = {
            "run_id": np.full(len(rates), data["run_id"][0], dtype=np.int32),
            "sec": np.full(len(rates), data["sec"][0], dtype=np.int64),
            "rates": np.array(rates, dtype=np.float64),
            "rates_err": np.array(rates_err, dtype=np.float64),
            "wp_tagging_efficiency": np.full(len(rates), wp_tagging_efficiency, dtype=np.float64)
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
        #     r"$\lambda = %.2f \pm %.2f~\mathrm{cps}$"
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

def gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu)/sigma)**2)

def plot_rate_hist(ax, data, color, title, interval=None):
    nbins = 50
    if interval is not None:
        bins = np.linspace(interval[0], interval[1], nbins + 1)
    else:
        bins = np.linspace(np.min(data), np.max(data), nbins + 1)
    h, bins, _ = ax.hist(data, bins=bins, alpha=0.6, color=color)
    e = np.sqrt(h)

    centers = 0.5 * (bins[1:] + bins[:-1])

    mask = h > 0
    x_fit = centers[mask]
    y_fit = h[mask]
    yerr_fit = e[mask]

    A0 = max(y_fit)
    mu0 = np.mean(data)
    sigma0 = np.std(data)
    p0 = [A0, mu0, sigma0]

    popt, pcov = curve_fit(gauss, x_fit, y_fit, p0=p0, sigma=yerr_fit, absolute_sigma=True)
    A, mu, sigma = popt
    A_err, mu_err, sigma_err = np.sqrt(np.diag(pcov))

    y_model = gauss(x_fit, *popt)
    chisq = np.sum(((y_model - y_fit) / yerr_fit)**2)
    ndf = len(y_fit) - len(popt)
    prob = chi2.sf(chisq, ndf)

    x_smooth = np.linspace(bins[0], bins[-1], 500)
    y_smooth = gauss(x_smooth, *popt)
    
    ax.plot(x_smooth, y_smooth, linestyle="--", linewidth=1.6, color="black", zorder=4)

    ax.set_xlabel(r"Muon rate (cps)")
    ax.set_ylabel(r"Entries")
    ax.set_title(title)

    text = (
        r"$P(\chi^2/\mathrm{ndf} = %.1f / %d) = %.3f$" "\n"
        r"$A = %.3f \pm %.3f$" "\n"
        r"$\mu = %.3f \pm %.3f$" "\n"
        r"$\sigma = %.3f \pm %.3f$"
    ) % (chisq, ndf, prob, A, A_err, mu, mu_err, sigma, sigma_err)
    ax.text(0.95, 0.95, text, transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')

def date_formatter(x, pos):
    date_str = mdates.num2date(x).strftime("%Y-%m-%d")
    return date_str.replace("-", r"\mbox{-}")

if __name__ == "__main__":
    args = parse_args()
    set_latex_style()

    if args.run > 0:
        calculate_muon_rate(args.run, plot=False)
        exit(0)
    
    if args.input == "":
        print("No input file specified")
        exit(1)

    run_ids = []
    timestamps = []
    rates_cd_wp = []
    rates_cd_only = []
    rates_wp_only = []
    rates_err_cd_wp = []
    rates_err_cd_only = []
    rates_err_wp_only = []

    file = uproot.open(args.input)
    tree = file["rates"]
    branches = ["run_id", "sec", "rates", "rates_err"]
    data = tree.arrays(branches, library="np")
    run_id = data["run_id"]
    sec = data["sec"]
    rates = data["rates"]
    rates_err = data["rates_err"]
    for i in range(len(run_id) // 3):
        run_ids.append(run_id[i * 3])
        timestamps.append(sec[i * 3])
        rates_cd_wp.append(rates[i * 3])
        rates_cd_only.append(rates[i * 3 + 1])
        rates_wp_only.append(rates[i * 3 + 2])
        rates_err_cd_wp.append(rates_err[i * 3])
        rates_err_cd_only.append(rates_err[i * 3 + 1])
        rates_err_wp_only.append(rates_err[i * 3 + 2])

    run_ids = np.array(run_ids)
    timestamps = np.array(timestamps)
    rates_cd_wp = np.array(rates_cd_wp)
    rates_cd_only = np.array(rates_cd_only)
    rates_wp_only = np.array(rates_wp_only)
    rates_err_cd_wp = np.array(rates_err_cd_wp)
    rates_err_cd_only = np.array(rates_err_cd_only)
    rates_err_wp_only = np.array(rates_err_wp_only)

    rates_total = rates_cd_wp + rates_cd_only + rates_wp_only
    err_total = np.sqrt(rates_err_cd_wp**2 + rates_err_cd_only**2 + rates_err_wp_only**2)

    mean_cd = np.mean(rates_cd_only)
    std_cd = np.std(rates_cd_only)

    exp = int(np.floor(np.log10(mean_cd)))
    mantissa_mean = mean_cd / 10**exp
    mantissa_std = std_cd / 10**exp

    cd_label = (
        rf"$\mathrm{{CD\ only}}: "
        rf"({mantissa_mean:.2f} \pm {mantissa_std:.2f})"
        rf"\times 10^{{{exp}}}\ \mathrm{{cps}}$"
    )

    config = [
        (rates_total, err_total, "#000000", rf"$\mathrm{{Total}}: {np.mean(rates_total):.2f} \pm {np.std(rates_total):.2f}\ \mathrm{{cps}}$"),
        (rates_cd_wp, rates_err_cd_wp, "#ffa500", rf"$\mathrm{{CD+WP}}: {np.mean(rates_cd_wp):.2f} \pm {np.std(rates_cd_wp):.2f}\ \mathrm{{cps}}$"),
        (rates_cd_only, rates_err_cd_only, "#1ea50d", cd_label),
        (rates_wp_only, rates_err_wp_only, "#3d80e6", rf"$\mathrm{{WP\ only}}: {np.mean(rates_wp_only):.2f} \pm {np.std(rates_wp_only):.2f}\ \mathrm{{cps}}$"),
]

    fig, ax = plt.subplots(figsize=(16, 6))

    for data, err, color, label in config:
        ax.errorbar(run_ids, data, yerr=err, fmt="o", color=color, markersize=5.0, capsize=0, elinewidth=1.0, markeredgecolor="k", markeredgewidth=0.5, label=label, zorder=3)
        mean_val = np.mean(data)
        ax.axhline(mean_val, color=color, linestyle="--", linewidth=2.0, zorder=2)

    ax.set_xlabel(r"Run Number")
    ax.set_ylabel(r"Muon rate (cps)")
    ax.set_ylim(-0.1, 10)

    ax.tick_params(direction='in', which='both', top=True, right=True)
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    # title_str = f"Run range: {min(run_ids)} - {max(run_ids)}" # exposure: {exposure_days:.1f} days"
    # ax.set_title(title_str, loc='right', color='grey', pad=20)

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, frameon=False, handletextpad=0.1)

    fig.tight_layout()
    fig.show()
    fig.savefig("muon_rate_vs_run.pdf")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    mean_cdwp = np.mean(rates_cd_wp)
    plot_rate_hist(axes[0], rates_cd_wp, "#ffa500", "CD+WP rate", interval=(mean_cdwp - 0.25, mean_cdwp + 0.25))

    mean_cd = np.mean(rates_cd_only)
    plot_rate_hist(axes[1], rates_cd_only, "#1ea50d", "CD Only rate", interval=(mean_cd - 0.005, mean_cd + 0.005))

    mean_wp = np.mean(rates_wp_only)
    plot_rate_hist(axes[2], rates_wp_only, "#3d80e6", "WP Only rate", interval=(mean_wp - 0.5, mean_wp + 0.5))

    fig.tight_layout()
    fig.show()

    dates = [datetime.fromtimestamp(ts, tz=timezone.utc) for ts in timestamps]

    fig, ax = plt.subplots(figsize=(16, 6))

    for data, err, color, label in config:
        ax.errorbar(
            dates,
            data,
            yerr=err,
            fmt="o",
            color=color,
            markersize=5.0,
            capsize=0,
            elinewidth=1.0,
            markeredgecolor="k",
            markeredgewidth=0.5,
            label=label,
            zorder=3
        )
        mean_val = np.mean(data)
        ax.axhline(mean_val, color=color, linestyle="--", linewidth=2.0, zorder=2)

    ax.xaxis.set_major_formatter(FuncFormatter(date_formatter))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())

    ax.set_ylabel(r"Muon rate (cps)")
    ax.set_ylim(-0.1, 10)

    ax.tick_params(direction='in', which='both', top=True, right=True)
    ax.minorticks_on()
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, frameon=False, handletextpad=0.1)

    # fig.autofmt_xdate()
    fig.tight_layout()
    fig.show()
    fig.savefig("muon_rate_vs_date.pdf")

    plt.show()