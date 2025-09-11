import uproot
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timezone

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

file1 = "~/Documents/test/crosscheck_sep8_vanessa.root"
file2 = "~/Documents/test/RUN.9737-9820.20250826-20250901.output.root"

tree1_name = "events"
tree2_name = "FirstCrossCheckAnalysis"

mapping = {
    "energy_p": ("e_p", np.arange(0.0, 12.0, 0.5), r"$E_{p}$ (MeV)"),
    "energy_d": ("e_d", np.arange(1.5, 3.0, 0.05), r"$E_{d}$ (MeV)"),
    "n_pe_p": ("totq_p", np.arange(0.0, 20000.0, 1000.0), r"$PEs$"),
    "n_pe_d": ("totq_d", np.arange(3500.0, 6500.0, 100.0), r"$PEs$"),
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
dr2 = np.linalg.norm(pos_p - pos_d, axis=1) / 1000.0

dt2 = []
for sp, np_, sd, nd in zip(data2["sec_p"], data2["nsec_p"], data2["sec_d"], data2["nsec_d"]):
    tp = TimeStamp(sp, np_)
    td = TimeStamp(sd, nd)
    dt2.append((td.to_sec() - tp.to_sec()) * 1000.0)
dt2 = np.array(dt2)

extra_vars = {
    "dr": (dr2, np.arange(0.0, 1.5, 0.05), r"$dr_{p2d}$ (m)"),
    "dt": (dt2, np.arange(0.0, 2.0, 0.1), r"$dt_{p2d}$ (ms)"),
}

# ---------------- Plotting ----------------
all_vars = list(mapping.keys()) + list(extra_vars.keys())
nvars = len(all_vars)

fig, axes = plt.subplots(2, nvars, figsize=(5*nvars, 10), sharex='col')

for i, var in enumerate(all_vars):
    if var in mapping:
        var2, bins, xlabel = mapping[var]
        vals1 = data1[var]
        vals2 = data2[var2]
    else:
        vals1 = data1[var]
        vals2, bins, xlabel = extra_vars[var]

    # --- Top: overlay histograms ---
    ax = axes[0, i]
    ax.hist(vals1, bins=bins, histtype="step", color="blue", label="Vanessa")
    ax.hist(vals2, bins=bins, histtype="step", color="red", label="Thomas")
    ax.set_title(f"{var} comparison")
    ax.set_ylabel("Entries")
    ax.legend()

    # --- Bottom: residuals (Vanessa - Thomas) ---
    counts1, _ = np.histogram(vals1, bins=bins)
    counts2, _ = np.histogram(vals2, bins=bins)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])

    ax_diff = axes[1, i]
    ax_diff.step(bin_centers, counts1 - counts2, where="mid", color="black")
    ax_diff.axhline(0, color="gray", linestyle="--")
    ax_diff.set_xlabel(xlabel)
    ax_diff.set_ylabel("Δ (V - T)")

plt.tight_layout()
plt.show()