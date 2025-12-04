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

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="Input filepath")
args = parser.parse_args()

with uproot.open(args.input) as f:
    tree_neu_all_oec = f["NeutronVetoStudy_All__Oec"]
    tree_cosmo_all_oec = f["CdWpCosmoStudy_All__Oec"]
    neu_all_oec = tree_neu_all_oec.arrays(library="np")
    cosmo_all_oec = tree_cosmo_all_oec.arrays(library="np")

fig_e_e, ax_e_e = plt.subplots(1, 1, figsize=(12, 10))
plt.hist(neu_all_oec["e_e"], bins=50)
set_plot_style(ax_e_e)
ax_e_e.set_title(r"$\text{Neutron energy}$")
ax_e_e.set_xlabel(r"$E_n$ (MeV)")
ax_e_e.set_ylabel(r"Entries")
fig_e_e.show()

d_np = np.sqrt(
    (neu_all_oec["posx_e"] - neu_all_oec["posx_p"])**2 +
    (neu_all_oec["posy_e"] - neu_all_oec["posy_p"])**2 +
    (neu_all_oec["posz_e"] - neu_all_oec["posz_p"])**2
) / 1000.0

fig_d_np, ax_d_np = plt.subplots(1, 1, figsize=(12, 10))
plt.hist(d_np, bins=50)
set_plot_style(ax_d_np)
ax_d_np.set_title(r"$\text{Distance neutron-prompt}$")
ax_d_np.set_xlabel(r"$d_{p-n}$ (m)")
ax_d_np.set_ylabel(r"Entries")
fig_d_np.show()

fig_e_e_d_np, ax_e_e_d_np = plt.subplots(1, 1, figsize=(12, 10))
h_e_e_d_np = plt.hist2d(neu_all_oec["e_e"], d_np, bins=(50, 50), cmin=1)
set_plot_style(ax_e_e_d_np)
ax_e_e_d_np.set_title(r"$\text{Neutron energy - Distance neutron-prompt}$")
ax_e_e_d_np.set_xlabel(r"$E_n$ (MeV)")
ax_e_e_d_np.set_ylabel(r"$d_{p-n}$ (m)")
cbar_e_e_d_np = fig_e_e_d_np.colorbar(h_e_e_d_np[3], ax=ax_e_e_d_np, label=r"Entries")
set_plot_colorbar_style(cbar_e_e_d_np.ax)
fig_e_e_d_np.show()

fig_e_e_e_p, ax_e_e_e_p = plt.subplots(1, 1, figsize=(12, 10))
h_e_e_e_p = plt.hist2d(neu_all_oec["e_e"], neu_all_oec["e_p"], bins=(50, 50), cmin=1)
set_plot_style(ax_e_e_e_p)
ax_e_e_e_p.set_title(r"$\text{Neutron energy - Prompt energy}$")
ax_e_e_e_p.set_xlabel(r"$E_n$ (MeV)")
ax_e_e_e_p.set_ylabel(r"$E_p$ (MeV)")
cbar_e_e_e_p = fig_e_e_e_p.colorbar(h_e_e_e_p[3], ax=ax_e_e_e_p, label=r"Entries")
set_plot_colorbar_style(cbar_e_e_e_p.ax)
fig_e_e_e_p.show()

dt_mu2p = [
    TimeStamp(sec, nsec).to_sec() for sec, nsec in zip(cosmo_all_oec["dt_mu2p_sec"], cosmo_all_oec["dt_mu2p_nsec"])
]

fig_dlat_p_dt_mu2p, ax_dlat_p_dt_mu2p = plt.subplots(1, 1, figsize=(12, 10))
h_dlat_p_dt_mu2p = plt.hist2d(cosmo_all_oec["dlat_p"] / 1000.0, dt_mu2p, bins=(50, 50), cmin=1)
set_plot_style(ax_dlat_p_dt_mu2p)
ax_dlat_p_dt_mu2p.set_title(r"$\text{Distance muon-prompt - Time muon-prompt}$")
ax_dlat_p_dt_mu2p.set_xlabel(r"$d_{\mu-p}$ (m)")
ax_dlat_p_dt_mu2p.set_ylabel(r"$t_{\mu-p}$ (s)")
cbar_dlat_p_dt_mu2p = fig_dlat_p_dt_mu2p.colorbar(h_dlat_p_dt_mu2p[3], ax=ax_dlat_p_dt_mu2p, label=r"Entries")
set_plot_colorbar_style(cbar_dlat_p_dt_mu2p.ax)
fig_dlat_p_dt_mu2p.show()

plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()

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
