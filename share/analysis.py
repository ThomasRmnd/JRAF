import argparse
from datetime import datetime, timezone

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator, AutoMinorLocator, MultipleLocator
import numpy as np
from scipy.optimize import curve_fit
import uproot

class TimeStamp:
    __slots__ = ('sec', 'nsec')

    def __init__(self, sec: int, nsec: int):
        self.sec = int(sec)
        self.nsec = int(nsec)
        self.normalize()

    def normalize(self):
        while (self.nsec < 0):
            self.nsec += 1_000_000_000
            self.sec -= 1
        while (self.nsec >= 1_000_000_000):
            self.nsec -= 1_000_000_000
            self.sec += 1

    def to_sec(self):
        return self.sec + self.nsec * 1e-9

    def to_nsec(self):
        return self.sec * 1_000_000_000 + self.nsec

    def __add__(self, other):
        return TimeStamp(self.sec + other.sec, self.nsec + other.nsec)

    def __sub__(self, other):
        return TimeStamp(self.sec - other.sec, self.nsec - other.nsec)

    def __lt__(self, other):
        return (self.sec, self.nsec) < (other.sec, other.nsec)

    def __le__(self, other):
        return (self.sec, self.nsec) <= (other.sec, other.nsec)

    def __gt__(self, other):
        return (self.sec, self.nsec) > (other.sec, other.nsec)

    def __ge__(self, other):
        return (self.sec, self.nsec) >= (other.sec, other.nsec)

    def __eq__(self, other):
        return (self.sec, self.nsec) == (other.sec, other.nsec)

    def __repr__(self):
        dt = datetime.fromtimestamp(self.sec, tz=timezone.utc)
        return dt.strftime(f"%Y-%m-%d %H:%M:%S.{self.nsec:09d} UTC")
    
mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "Times"],
    # small preamble to support \text and amsmath usage
    "text.latex.preamble": r"\usepackage{amsmath}\usepackage{siunitx}",
    "font.size": 14,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "figure.dpi": 100,
})
    
def set_plot_style(ax):
    """
    Apply uniform 'publication-like' style:
      - LaTeX font
      - Major ticks adapting to axis range
      - Minor ticks automatically subdivided
      - Inward ticks on all sides
    """
    # Major ticks: auto locator
    ax.xaxis.set_major_locator(AutoLocator())
    ax.yaxis.set_major_locator(AutoLocator())

    # Minor ticks: automatically subdivide (default: 4 per major interval)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    # Tick params: inward, both sides
    ax.tick_params(which="both", direction="in", top=True, right=True,
                   length=4, width=1, labelsize=12)
    ax.tick_params(which="major", length=6, width=1.2)
    ax.tick_params(which="minor", length=3, width=1)

def set_plot_colorbar_style(ax):
    """
    Apply uniform 'publication-like' style:
      - LaTeX font
      - Major ticks adapting to axis range
      - Minor ticks automatically subdivided
      - Inward ticks on all sides
    """
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    ax.tick_params(axis="y", which="both", direction="in", right=True, left=True)
    
def analyze_daqtree(filename):
    """Special analysis for Thomas' DAQTree"""
    with uproot.open(filename) as f:
        if "DAQTree" not in f:
            print("[INFO] No DAQTree found in file")
            return

        tree = f["DAQTree"]

        daq_sec = tree["daq_sec"].array(library="np")
        daq_nsec = tree["daq_nsec"].array(library="np")
        muveto_sec = tree["muveto_sec"].array(library="np")
        muveto_nsec = tree["muveto_nsec"].array(library="np")

        daq_ts = daq_sec.astype(np.float64) + daq_nsec * 1e-9
        muveto_ts = muveto_sec.astype(np.float64) + muveto_nsec * 1e-9

        daq_sum = daq_ts.sum()
        muveto_sum = muveto_ts.sum()

        daq_hours = daq_sum / 3600
        daq_days = daq_sum / 86400
        muveto_hours = muveto_sum / 3600
        muveto_days = muveto_sum / 86400

        perc = 100.0 * muveto_sum / daq_sum if daq_sum > 0 else float("nan")

        print("=== DAQTree Analysis ===")
        print(f"DAQ total time     : {daq_sum:.3e} s = {daq_hours:.3f} h = {daq_days:.3f} d")
        print(f"MuVeto total time  : {muveto_sum:.3e} s = {muveto_hours:.3f} h = {muveto_days:.3f} d")
        print(f"MuVeto / DAQ ratio : {perc:.3f} %")
        print("========================")
    
class CodeAnalysis:
    def __init__(self, name: str, config: dict):
        self.name = name
        self.filename = config["filename"]
        self.treename = config["treename"]
        self.varmap = config.get("varmap", {})
        self.derived = config.get("derived", {})
        self.data = {}

    def load(self):
        with uproot.open(self.filename) as f:
            tree = f[self.treename]
            for var, treename in self.varmap.items():
                self.data[var] = tree[treename].array(library="np")

    def get(self, var):
        if var in self.data:
            return self.data[var]

        if var in self.derived:
            self.data[var] = self.derived[var](self)
            return self.data[var]

        raise KeyError(f"{var} not available for {self.name}")

def dr_thomas(ana):
    dx = ana.get("posx_p") - ana.get("posx_d")
    dy = ana.get("posy_p") - ana.get("posy_d")
    dz = ana.get("posz_p") - ana.get("posz_d")
    return np.sqrt(dx**2 + dy**2 + dz**2)

def dt_thomas(ana):
    sec_p, nsec_p = ana.get("sec_p"), ana.get("nsec_p")
    sec_d, nsec_d = ana.get("sec_d"), ana.get("nsec_d")

    ts_p = np.array([TimeStamp(s, ns) for s, ns in zip(sec_p, nsec_p)])
    ts_d = np.array([TimeStamp(s, ns) for s, ns in zip(sec_d, nsec_d)])

    dt_vals = np.array([(td - tp).to_sec() for tp, td in zip(ts_p, ts_d)])
    return dt_vals


def rho_p_thomas(ana):
    return np.sqrt(ana.get("posx_p")**2 + ana.get("posy_p")**2)

def z_p_thomas(ana):
    return ana.get("posz_p")

def rho_d_thomas(ana):
    return np.sqrt(ana.get("posx_d")**2 + ana.get("posy_d")**2)

def z_d_thomas(ana):
    return ana.get("posz_d")
    
parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="Input filepath")
parser.add_argument("--input-vanessa", type=str, default="", help="Input Vanessa filepath")
parser.add_argument("--input-cristobal", type=str, default="", help="Input Cristobal filepath")
args = parser.parse_args()

analyzer_configs = {
    "Thomas": {
        "filename": args.input,
        "treename": "FirstCrossCheckAnalysis",
        "varmap": {
            "e_p": "e_p",
            "e_d": "e_d",
            "totq_p": "totq_p",
            "totq_d": "totq_d",
            "posx_p": "posx_p",
            "posy_p": "posy_p",
            "posz_p": "posz_p",
            "posx_d": "posx_d",
            "posy_d": "posy_d",
            "posz_d": "posz_d",
            "sec_p": "sec_p",
            "nsec_p": "nsec_p",
            "sec_d": "sec_d",
            "nsec_d": "nsec_d",
        },
        "derived": {
            "dr": dr_thomas,
            "dt": dt_thomas,
            "rho_p": rho_p_thomas,
            "rho_d": rho_d_thomas,
            "z_p": z_p_thomas,
            "z_d": z_d_thomas
        }
    },
    "Vanessa": {
        "filename": args.input_vanessa,
        "treename": "events",
        "varmap": {
            "e_p": "energy_p",
            "e_d": "energy_d",
            "totq_p": "n_pe_p",
            "totq_d": "n_pe_d",
        }
    },
    "Cristobal": {
        "filename": args.input_cristobal,
        "treename": "ibds",
        "varmap": {
            "e_p": "p_energy",
            "e_d": "d_energy",
            "totq_p": "p_charge",
            "totq_d": "d_charge",
        }
    }
}

variables = {
    "e_p": {
        "bins": np.linspace(0.0, 12.0, 51), 
        "title": r"$\text{Prompt energy}$", 
        "xlabel": r"$E_p$ (MeV)"
    },
    "e_d": {
        "bins": np.linspace(1.5 ,3.0 , 51), 
        "title": r"$\text{Delayed energy}$", 
        "xlabel": r"$E_d$ (MeV)"
    },
    "totq_p": {
        "bins": np.linspace(500.0, 21000.0, 51), 
        "title": r"$\text{Prompt total charge}$", 
        "xlabel": r"$\text{Prompt PEs}$"
    },
    "totq_d": {
        "bins": np.linspace(3500.0, 6500.0, 51), 
        "title": r"$\text{Delayed total charge}$", 
        "xlabel": r"$\text{Delayed PEs}$"
    },
    "dr": {
        "bins": np.linspace(0.0, 1.5, 51),
        "title": r"$\text{Prompt-delayed distance}$",
        "xlabel": r"$\Delta r_{p-d}$ (m)"
    },
    "dt": {
        "bins": np.linspace(0.0, 2.0, 51),
        "title": r"$\text{Prompt-delayed time}$",
        "xlabel": r"$\Delta t_{p-d}$ (ms)"
    },
}

analyzers = [
    CodeAnalysis(name, config) for name, config in analyzer_configs.items()
    if config["filename"] != ""
]
for ana in analyzers:
    ana.load()

# if args.input:
#     analyze_daqtree(args.input)

# ---------------- Fitting function ----------------
def exp_decay(x, A, tau):
    return A * np.exp(-x / tau)

def plot_group_comparison(analyzers, var_list):
    colors = [
        "tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", 
        "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"
    ]
    for var in var_list:
        if var in variables:
            bins = variables[var]["bins"]
            title = variables[var]["title"]
            xlabel = variables[var]["xlabel"]

        ref_ana = analyzers[0]
        vals_ref = ref_ana.get(var)
        counts_ref, _ = np.histogram(vals_ref, bins=bins)
        bin_centers = 0.5 * (bins[1:] + bins[:-1])

        fig, (ax_top, ax_diff) = plt.subplots(
            2, 1, figsize=(8, 8), sharex=True,
            gridspec_kw={'height_ratios': [3, 1]}
        )

        for ana, color in zip(analyzers, colors):
            vals = ana.get(var)
            ax_top.hist(vals, bins=bins, histtype="step", color=color, linewidth=1.6, label=ana.name)

            counts, _ = np.histogram(vals, bins=bins)
            if ana is not ref_ana:
                ax_diff.step(bin_centers, counts - counts_ref,
                             where="mid", color=color, linewidth=1.6, label=f"{ana.name}-{ref_ana.name}")

        set_plot_style(ax_top)
        set_plot_style(ax_diff)

        # Labels (use LaTeX expressions)
        ax_top.set_title(title)
        ax_top.set_ylabel(r"Entries")
        ax_top.legend(frameon=False)

        ax_diff.axhline(0, color="gray", linestyle="--")
        ax_diff.set_xlabel(xlabel)
        ax_diff.set_ylabel(r"$\Delta$ Entries")

        plt.tight_layout(pad=0.6)

def plot_geometry_and_time(analyzers):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    ana = analyzers[0]  # e.g. Thomas

    # --- rho2_p vs z_p ---
    ax = axes[0, 0]
    set_plot_style(ax)
    rho2 = ana.get("rho_p")**2 / 1e6 # m^2
    z = ana.get("posz_p") / 1000.0  # m
    h = ax.hist2d(rho2, z, bins=(50, 50), cmap="viridis", cmin=1)
    ax.set_xlabel(r"$\rho^2_p$ (m$^2$)")
    ax.set_ylabel(r"$z_p$ (m)")
    cbar = fig.colorbar(h[3], ax=ax, label=r"Entries")
    set_plot_colorbar_style(cbar.ax)   # make colorbar ticks consistent

    # --- rho2_d vs z_d ---
    ax = axes[0, 1]
    set_plot_style(ax)
    rho2 = ana.get("rho_d")**2 / 1e6 # m^2
    z = ana.get("posz_d") / 1000.0  # m
    h = ax.hist2d(rho2, z, bins=(50, 50), cmap="viridis", cmin=1)
    ax.set_xlabel(r"$\rho^2_d$ (m$^2$)")
    ax.set_ylabel(r"$z_d$ (m)")
    cbar = fig.colorbar(h[3], ax=ax, label=r"Entries")
    set_plot_colorbar_style(cbar.ax)

    # --- dr histogram ---
    ax = axes[1, 0]
    set_plot_style(ax)
    dr_vals = ana.get("dr") / 1000.0 # m
    ax.hist(dr_vals, bins=np.linspace(0, 1.5, 51), histtype="step", color="tab:blue")
    ax.set_xlabel(r"$\Delta r_{p-d}$ (m)")
    ax.set_ylabel(r"Entries")

    # --- dt histogram with exponential fit ---
    ax = axes[1, 1]
    set_plot_style(ax)
    dt_vals = ana.get("dt") * 1000.0 # ms
    counts, edges = np.histogram(dt_vals, bins=np.linspace(0.0, 2.0, 51))
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    ax.errorbar(bin_centers, counts, yerr=np.sqrt(counts), fmt='o', color="tab:blue")

    mask = counts > 0
    if mask.sum() > 2:
        popt, pcov = curve_fit(exp_decay, bin_centers[mask], counts[mask], p0=(counts.max(), 0.5))
        A_fit, tau_fit_ms = popt
        tau_err_ms = np.sqrt(np.diag(pcov))[1]
        tau_fit_us = tau_fit_ms * 1000.0
        tau_err_us = tau_err_ms * 1000.0
        ax.plot(bin_centers, exp_decay(bin_centers, *popt), "b--",
            label=fr"$\tau = {tau_fit_us:.1f} \pm {tau_err_us:.1f}\,\mu\mathrm{{s}}$")
        ax.legend(frameon=False)

    ax.set_xlabel(r"$\Delta t_{p-d}$ (ms)")
    ax.set_ylabel(r"Entries")

    plt.tight_layout(pad=0.6)

plot_group_comparison(analyzers, ["e_p", "e_d", "totq_p", "totq_d"])
plot_geometry_and_time(analyzers)

plt.show()