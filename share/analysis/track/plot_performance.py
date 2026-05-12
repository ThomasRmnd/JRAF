import argparse

import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FuncFormatter
from matplotlib.transforms import offset_copy
import numpy as np
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

set_latex_style()

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, required=True, help="Input filepath")
parser.add_argument("--label", type=str, required=True, help="Label name")
args = parser.parse_args()

with uproot.open(args.input) as f:
    t              = f["performance"]
    angle          = t["angle"]
    dist_mid_point = t["dist_mid_point"]
    dist_center    = t["dist_center"]

angle_xmin = 0.0
angle_xmax = 5.0
angle_nbin = 50

angle_full_xmin = 0.0
angle_full_xmax = 180.0
angle_full_nbin = int(angle_nbin * (angle_full_xmax - angle_full_xmin) / (angle_xmax - angle_xmin))

angle_color          = "#648fff"
dist_mid_point_color = "#ff6464"

angle_label          = r"$\alpha$ (deg)"
dist_mid_point_label = r"$d_{\mathrm{mid}}$ (m)"
dist_center_label    = r"$R_{\mu}^2$ (m$^2$)"

dist_mid_point_xmin = 0
dist_mid_point_xmax = 2.0
dist_mid_point_nbin = 50

dist_mid_point_full_xmin = 0.0
dist_mid_point_full_xmax = 40.0
dist_mid_point_full_nbin = int(dist_mid_point_nbin * (dist_mid_point_full_xmax - dist_mid_point_full_xmin) / (dist_mid_point_xmax - dist_mid_point_xmin))

fig, ax = plt.subplots(figsize=(7, 6))

bins = np.linspace(angle_full_xmin, angle_full_xmax, angle_full_nbin + 1)

hist, edges = np.histogram(angle, bins=bins, density=True)
p68 = np.quantile(angle, 0.68)
ax.fill_between(bins, np.r_[hist, hist[-1]], step="post", color=angle_color, alpha=0.075, zorder=1)
ax.step(bins, np.r_[hist, hist[-1]], where="post", color=angle_color, linestyle="-", linewidth=1.5, zorder=2, label=rf"{args.input}" + r", $68\%$ = " + f"{p68:.1f} deg")

ax.set_xlabel(angle_label)
ax.set_ylabel(r"P.D.F.")
ax.set_xlim(angle_xmin, angle_xmax)
ax.set_ylim(bottom=0.0)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
ax.legend(loc="upper right")

fig.tight_layout()
fig.show()

fig, ax = plt.subplots(figsize=(7, 6))

bins = np.linspace(dist_mid_point_full_xmin, dist_mid_point_full_xmax, dist_mid_point_full_nbin + 1)

hist, edges = np.histogram(dist_mid_point, bins=bins, density=True)
p68 = np.quantile(dist_mid_point, 0.68)
ax.fill_between(bins, np.r_[hist, hist[-1]], step="post", color=dist_mid_point_color, alpha=0.075, zorder=1)
ax.step(bins, np.r_[hist, hist[-1]], where="post", color=dist_mid_point_color, linestyle="-", linewidth=1.5, zorder=2, label=rf"{args.input}" + r", $68\%$ = " + f"{p68:.2f} m")

ax.set_xlabel(dist_mid_point_label)
ax.set_ylabel(r"P.D.F.")
ax.set_xlim(dist_mid_point_xmin, dist_mid_point_xmax)
ax.set_ylim(bottom=0.0)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
ax.legend(loc="upper right")

fig.tight_layout()
fig.show()