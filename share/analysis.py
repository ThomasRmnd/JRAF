import uproot
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timezone
from scipy.optimize import curve_fit

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

# Input ROOT files
file1 = "~/Documents/test/crosscheck_sep8_vanessa.root"
file2 = "~/Documents/test/RUN.20250826_20250907.summary.root"

tree1_name = "events"
tree2_name = "FirstCrossCheckAnalysis"

mapping = {
    "energy_p": ("e_p", np.linspace(0.0, 12.0, 51), r"$E_{p}$ (MeV)"),
    "energy_d": ("e_d", np.linspace(1.5, 3.0, 51), r"$E_{d}$ (MeV)"),
    "n_pe_p": ("totq_p", np.linspace(500.0, 21000.0, 51), "Prompt PEs"),
    "n_pe_d": ("totq_d", np.linspace(3500.0, 6500.0, 51), "Delayed PEs"),
}

# ---------------- Load ROOT files ----------------
with uproot.open(file1) as f1, uproot.open(file2) as f2:
    tree1 = f1[tree1_name]
    tree2 = f2[tree2_name]

    data1 = tree1.arrays(list(mapping.keys()) + ["dr", "dt"], library="np")
    data2 = tree2.arrays([
        *[v[0] for v in mapping.values()],
        "posx_p","posy_p","posz_p",
        "posx_d","posy_d","posz_d",
        "sec_p","nsec_p","sec_d","nsec_d"
    ], library="np")

# ---------------- Derived variables ----------------
pos_p = np.vstack([data2["posx_p"], data2["posy_p"], data2["posz_p"]]).T
pos_d = np.vstack([data2["posx_d"], data2["posy_d"], data2["posz_d"]]).T
dr2 = np.linalg.norm(pos_p - pos_d, axis=1) / 1000.0  # m

dt2 = []
for sp, np_, sd, nd in zip(data2["sec_p"], data2["nsec_p"], data2["sec_d"], data2["nsec_d"]):
    tp = TimeStamp(sp, np_)
    td = TimeStamp(sd, nd)
    dt2.append((td.to_sec() - tp.to_sec()) * 1000.0)  # ms
dt2 = np.array(dt2)

rho2_p = data2["posx_p"]**2 + data2["posy_p"]**2
rho2_d = data2["posx_d"]**2 + data2["posy_d"]**2
z_p = data2["posz_p"]
z_d = data2["posz_d"]

extra_vars = {
    "dr": (dr2, np.linspace(0.0, 1.5, 51), r"$\Delta r_{p-d}$ (m)"),
    "dt": (dt2, np.linspace(0.0, 2.0, 51), r"$\Delta t_{p-d}$ (ms)"),
}

ex_extra_vars = {
    "rho2_p": (rho2_p, np.linspace(0, 325e6, 51), r"Prompt $\rho^2$ (mm)"),
    "rho2_d": (rho2_d, np.linspace(0, 325e6, 51), r"Delayed $\rho^2$ (mm)"),
    "z_p": (z_p, np.linspace(-20e3, 20e3, 51), r"Prompt $z$ (mm)"),
    "z_d": (z_d, np.linspace(-20e3, 20e3, 51), r"Delayed $z$ (mm)"),
}

print("Number of events from Vanessa: ", len(data1["energy_p"]))
print("Number of events from Thomas: ", len(data2["e_p"]))

# ---------------- Fitting function ----------------
def exp_decay(x, A, tau):
    return A * np.exp(-x / tau)

# ---------------- Comparison plots (file1 vs file2) ----------------
def plot_group_comparison(var_list, title):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex='col')
    axes = axes.reshape(2, 2)
    for i, var in enumerate(var_list):
        row, col = divmod(i, 2)
        if var in mapping:
            var2, bins, xlabel = mapping[var]
            vals1 = data1[var]
            vals2 = data2[var2]
        else:
            vals1 = data1[var]
            vals2, bins, xlabel = extra_vars[var]

        # Top: overlay
        ax = axes[0, col]
        ax.hist(vals1, bins=bins, histtype="step", color="blue", label="Vanessa")
        ax.hist(vals2, bins=bins, histtype="step", color="red", label="Thomas")
        ax.set_title(f"{var}")
        ax.set_ylabel("Entries")
        ax.legend()

        # Bottom: residual
        counts1, _ = np.histogram(vals1, bins=bins)
        counts2, _ = np.histogram(vals2, bins=bins)
        bin_centers = 0.5 * (bins[1:] + bins[:-1])

        ax_diff = axes[1, col]
        ax_diff.step(bin_centers, counts1 - counts2, where="mid", color="black")
        ax_diff.axhline(0, color="gray", linestyle="--")
        ax_diff.set_xlabel(xlabel)
        ax_diff.set_ylabel("Δ (V - T)")

    fig.suptitle(title, fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])


# ---------------- File2-only plots (one row only, with rho² vs z) ----------------
def plot_group_data2(var_list, title):
    plt.style.use("default")  # white background

    fig, axes = plt.subplots(1, len(var_list), figsize=(6 * len(var_list), 6))

    # If only one variable, axes is not a list
    if len(var_list) == 1:
        axes = [axes]

    for i, var in enumerate(var_list):
        ax = axes[i]

        if var in mapping:
            var2, bins, xlabel = mapping[var]
            vals2 = data2[var2]

            counts, bins, _ = ax.hist(vals2, bins=bins, histtype="step", color="tab:blue")

            # If dt: fit exponential (convert ms -> µs for tau)
            if var == "dt":
                bin_centers = 0.5 * (bins[:-1] + bins[1:])
                mask = (counts > 0)
                if mask.sum() > 2:
                    popt, _ = curve_fit(exp_decay, bin_centers[mask], counts[mask],
                                        p0=(counts.max(), 0.5))
                    A_fit, tau_fit_ms = popt
                    tau_fit_us = tau_fit_ms * 1000.0
                    ax.plot(bin_centers, exp_decay(bin_centers, *popt), "b--",
                            label=f"Fit: A exp(-x/τ)\nτ = {tau_fit_us:.1f} µs")

            ax.set_xlabel(xlabel)
            ax.set_ylabel("Entries")

        elif var in extra_vars:
            vals2, bins, xlabel = extra_vars[var]

            if var == "dt":
                # Histogram counts and bin centers
                counts, edges = np.histogram(vals2, bins=bins)
                bin_centers = 0.5 * (edges[:-1] + edges[1:])
                errors = np.sqrt(counts)  # Poisson errors

                # Plot with error bars
                ax.errorbar(bin_centers, counts, yerr=errors, fmt='o', color="tab:blue")

                # Fit exponential decay
                mask = counts > 0
                if mask.sum() > 2:
                    popt, pcov = curve_fit(exp_decay, bin_centers[mask], counts[mask], p0=(counts.max(), 0.5))
                    A_fit, tau_fit_ms = popt
                    tau_err_ms = np.sqrt(np.diag(pcov))[1]  # uncertainty of tau in ms
                    tau_fit_us = tau_fit_ms * 1000.0       # convert to µs
                    tau_err_us = tau_err_ms * 1000.0       # convert to µs
                    ax.plot(bin_centers, exp_decay(bin_centers, *popt), "b--",
                        label=f"Fit: A exp(-x/τ)\nτ = {tau_fit_us:.1f} ± {tau_err_us:.1f} µs")
                    ax.legend()
            else:
                # Other variables: plain histogram
                ax.hist(vals2, bins=bins, histtype="step", color="tab:blue")

            ax.set_xlabel(xlabel)
            ax.set_ylabel("Entries")

        elif var == "rho2_p" or var == "rho2_d":
            # Prompt
            rho2, rho_bins, rho_label = ex_extra_vars[var]
            rho2 /= 1e6  # m²
            rho_bins /= 1e6  # m²
            rho_label = rho_label.replace("mm", "m")
            if var == "rho2_p":
                z, z_bins, z_label = ex_extra_vars["z_p"]
            if var == "rho2_d":
                z, z_bins, z_label = ex_extra_vars["z_d"]
            z /= 1000.0  # m
            z_bins /= 1000.0  # m
            z_label = z_label.replace("mm", "m")

            # Make white background
            ax.set_facecolor("white")
            fig.patch.set_facecolor("white")

            # 2D histogram
            h = ax.hist2d(rho2, z, bins=(rho_bins, z_bins), cmap="viridis", cmin=1)

            ax.set_xlabel(rho_label)
            ax.set_ylabel(z_label)

            # Add colorbar
            cbar = plt.colorbar(h[3], ax=ax)
            cbar.set_label("Counts")

    plt.tight_layout(rect=[0, 0, 1, 0.95])


# ---------------- Make all figures ----------------
plot_group_comparison(["energy_p", "energy_d"], "Energy Comparison")
plot_group_comparison(["n_pe_p", "n_pe_d"], "PEs Comparison")
plot_group_comparison(["dr", "dt"], "Distance & Time Comparison")

plot_group_data2(["energy_p", "energy_d"], "Energy Distributions")
plot_group_data2(["n_pe_p", "n_pe_d"], "PEs Distributions")
plot_group_data2(["dr", "dt"], "Distance & Time")
plot_group_data2(["rho2_p"], "Geometry")
plot_group_data2(["rho2_d"], "Geometry")

plt.show()