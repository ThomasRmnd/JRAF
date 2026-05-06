import argparse
from datetime import datetime, timedelta, timezone

import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FuncFormatter
from matplotlib.transforms import offset_copy
import numdifftools as nd
import numpy as np
from scipy.linalg import inv
from scipy.optimize import minimize
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
    parser.add_argument("--data", type=str, help="Cosmo data filepath")
    parser.add_argument("--simulation", type=str, help="Cosmo simu filepath")
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

def build_chi2(signal_e_p, background_e_p, mc_e_p, mc_element, bins_e, signal_dt_mu2p, bins_dt, fixed = None):
    fixed = fixed or {}

    S_sig, _      = np.histogram(signal_e_p,     bins=bins_e)
    S_bkg, _      = np.histogram(background_e_p, bins=bins_e)

    li9_mask      = mc_element == "Li9"
    he8_mask      = mc_element == "He8"
    S_Li9, _      = np.histogram(mc_e_p[li9_mask], bins=bins_e)
    S_He8, _      = np.histogram(mc_e_p[he8_mask], bins=bins_e)

    def safe_norm(h):
        n = h.sum()
        return h / n if n > 0 else h.astype(float)

    S_bkg_n  = safe_norm(S_bkg.astype(float))
    S_Li9_n  = safe_norm(S_Li9.astype(float))
    S_He8_n  = safe_norm(S_He8.astype(float))

    N_sig = S_sig.sum()   # total counts in signal region — kept fixed

    dt_centres = 0.5 * (bins_dt[:-1] + bins_dt[1:])
    dt_width   = np.diff(bins_dt)
    dt_mask    = dt_centres > 0.060

    dt_centres = dt_centres[dt_mask]
    dt_width   = dt_width[dt_mask]

    dt_sig_full, _ = np.histogram(signal_dt_mu2p, bins=bins_dt)
    dt_sig         = dt_sig_full[dt_mask]

    PARAM_NAMES   = ["f_bkg", "f_Li9", "tau_Li9", "tau_He8"]
    DEFAULT_VALS  = { "f_bkg": 0.25, "f_Li9": 0.90,
                      "tau_Li9": 0.257, "tau_He8": 0.172 }
    BOUNDS_MAP    = { "f_bkg": (0, 1), "f_Li9": (0, 1),
                      "tau_Li9": (0.05, 1.0), "tau_He8": (0.05, 1.0) }

    free_names  = [p for p in PARAM_NAMES if p not in fixed]
    free_bounds = [BOUNDS_MAP[p] for p in free_names]
    free_x0     = [DEFAULT_VALS[p] for p in free_names]

    def unpack(free_vals):
        d = dict(fixed)
        d.update(dict(zip(free_names, free_vals)))
        return d["f_bkg"], d["f_Li9"], d["tau_Li9"], d["tau_He8"]

    def chi2(free_vals):
        f_bkg, f_Li9, tau_Li9, tau_He8 = unpack(free_vals)

        # if not (0 <= f_bkg  <= 1): return 1e10
        # if not (0 <= f_Li9  <= 1): return 1e10
        # if tau_Li9 <= 0 or tau_He8 <= 0: return 1e10

        E_model = N_sig * (f_bkg * S_bkg_n + (1 - f_bkg) * (f_Li9 * S_Li9_n + (1 - f_Li9) * S_He8_n))

        sigma_E = np.sqrt(f_bkg * S_bkg.astype(float) + (1 - f_bkg) * (f_Li9 * S_Li9.astype(float) + (1 - f_Li9) * S_He8.astype(float)))
        sigma_E = np.where(sigma_E > 0, sigma_E, 1.0)

        chi2_E = np.sum(((S_sig - E_model) / sigma_E) ** 2)

        exp_Li9 = np.exp(-np.abs(dt_centres) / tau_Li9) * dt_width
        exp_He8 = np.exp(-np.abs(dt_centres) / tau_He8) * dt_width
        flat    = np.ones_like(dt_centres)  * dt_width

        def safe_norm_arr(a):
            s = a.sum()
            return a / s if s > 0 else a

        dt_model = N_sig * (f_bkg * safe_norm_arr(flat) + (1 - f_bkg) * (f_Li9 * safe_norm_arr(exp_Li9) + (1 - f_Li9) * safe_norm_arr(exp_He8)))

        sigma_dt_shape = (f_bkg * safe_norm_arr(flat) * N_sig + (1 - f_bkg) * f_Li9 * safe_norm_arr(exp_Li9) * N_sig + (1 - f_bkg) * (1 - f_Li9) * safe_norm_arr(exp_He8) * N_sig)
        sigma_dt = np.where(sigma_dt_shape > 0, np.sqrt(sigma_dt_shape), 1.0)

        chi2_dt = np.sum(((dt_sig - dt_model) / sigma_dt) ** 2)

        return chi2_E + chi2_dt

    return chi2, free_names, free_bounds, free_x0, S_sig, S_bkg, S_Li9, S_He8, dt_sig

def run_fit(signal_e_p, signal_dt_mu2p, background_e_p, mc_e_p, mc_element, bins_e, background_dt_mu2p, bins_dt, fixed = None):
    fixed = fixed or {}

    chi2, free_names, free_bounds, free_x0, S_sig, S_bkg, S_Li9, S_He8, dt_sig = build_chi2(
        signal_e_p, background_e_p, 
        mc_e_p, mc_element, bins_e,
        signal_dt_mu2p, bins_dt,
        fixed=fixed
    )

    result = minimize(chi2, free_x0, method="L-BFGS-B", bounds=free_bounds,
                      options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 100_000})

    # ── Uncertainties on free parameters only ─────────────────────────────────
    H      = nd.Hessian(chi2)(result.x)
    cov    = 2.0 * inv(H)
    errors = np.sqrt(np.diag(cov))

    ndof = (len(S_sig) + len(dt_sig)) - len(free_names)

    print(f"Fit converged : {result.success}")
    print(f"chi2 / ndof   = {result.fun:.1f} / {ndof} = {result.fun/ndof:.3f}")

    FMT = {"f_bkg": (".4f", ""),       "f_Li9":   (".4f", ""),
           "tau_Li9": (".1f", " ms"),  "tau_He8": (".1f", " ms")}
    SCALE = {"f_bkg": 1, "f_Li9": 1, "tau_Li9": 1e3, "tau_He8": 1e3}
    PDG   = {"tau_Li9": "PDG: 257 ms", "tau_He8": "PDG: 172 ms"}

    for name, val, err in zip(free_names, result.x, errors):
        fmt, unit = FMT[name]
        sc        = SCALE[name]
        pdg_note  = f"  ({PDG[name]})" if name in PDG else ""
        print(f"  {name:8s} = {val*sc:{fmt}}{unit}  +/-  {err*sc:{fmt}}{unit}{pdg_note}")

    for name, val in fixed.items():
        fmt, unit = FMT[name]
        sc        = SCALE[name]
        print(f"  {name:8s} = {val*sc:{fmt}}{unit}  [fixed]")

    # Rebuild full result.x in canonical order for plot_fit
    PARAM_NAMES  = ["f_bkg", "f_Li9", "tau_Li9", "tau_He8"]
    DEFAULT_VALS = dict(fixed)
    DEFAULT_VALS.update(dict(zip(free_names, result.x)))
    full_params  = np.array([DEFAULT_VALS[p] for p in PARAM_NAMES])

    full_errors  = {}
    for name, err in zip(free_names, errors):
        full_errors[name] = err
    for name in fixed:
        full_errors[name] = 0.0
    full_errors_arr = np.array([full_errors[p] for p in PARAM_NAMES])

    result.x = full_params   # patch so plot_fit stays unchanged
    return result, full_errors_arr, S_sig, S_bkg, S_Li9, S_He8, dt_sig, chi2, free_names

def plot_fit(result, errors, S_sig, S_bkg, S_Li9, S_He8, dt_sig,
             bins_e, bins_dt, chi2_func, free_names, fixed=None):

    fixed = fixed or {}
    f_bkg, f_Li9, tau_Li9, tau_He8 = result.x
    N_sig = S_sig.sum()

    def safe_norm(h):
        n = h.sum()
        return h / n if n > 0 else h.astype(float)

    # ── Reconstructed energy components ───────────────────────────────────────
    E_centres = 0.5 * (bins_e[:-1] + bins_e[1:])

    comp_bkg_E = N_sig * f_bkg             * safe_norm(S_bkg.astype(float))
    comp_Li9_E = N_sig * (1 - f_bkg)       * f_Li9       * safe_norm(S_Li9.astype(float))
    comp_He8_E = N_sig * (1 - f_bkg)       * (1 - f_Li9) * safe_norm(S_He8.astype(float))
    total_E    = comp_bkg_E + comp_Li9_E + comp_He8_E

    # ── Reconstructed dt components ───────────────────────────────────────────
    dt_centres_full = 0.5 * (bins_dt[:-1] + bins_dt[1:])
    dt_width_full   = np.diff(bins_dt)
    dt_mask         = dt_centres_full > 0.060
    dt_centres      = dt_centres_full[dt_mask]
    dt_width        = dt_width_full[dt_mask]

    def safe_norm_arr(a):
        s = a.sum()
        return a / s if s > 0 else a

    exp_Li9 = np.exp(-dt_centres / tau_Li9) * dt_width
    exp_He8 = np.exp(-dt_centres / tau_He8) * dt_width
    flat    = np.ones_like(dt_centres) * dt_width

    comp_bkg_dt = N_sig * f_bkg       * safe_norm_arr(flat)
    comp_Li9_dt = N_sig * (1 - f_bkg) * f_Li9       * safe_norm_arr(exp_Li9)
    comp_He8_dt = N_sig * (1 - f_bkg) * (1 - f_Li9) * safe_norm_arr(exp_He8)
    total_dt    = comp_bkg_dt + comp_Li9_dt + comp_He8_dt

    # ── Shared style ──────────────────────────────────────────────────────────
    kw_bkg  = dict(color="#FFA500", lw=1.5)
    kw_Li9  = dict(color="#00BFFF", lw=1.5)
    kw_He8  = dict(color="#9B59B6", lw=1.5)
    kw_sum  = dict(color="black",   lw=1.5, linestyle="dotted")
    kw_data = dict(fmt="ko", ms=4, capsize=3, lw=1.0, zorder=5)

    def step_plot(ax, bins, values, **kwargs):
        ax.step(bins, np.r_[values, values[-1]], where="post", **kwargs)

    def fill_plot(ax, bins, values, color):
        ax.fill_between(bins, np.r_[values, values[-1]],
                        step="post", color=color, alpha=0.15)

    # ── Energy plot ───────────────────────────────────────────────────────────
    fig_e, ax_e = plt.subplots(figsize=(7, 6))
    ax_e.errorbar(E_centres, S_sig, yerr=np.sqrt(S_sig.astype(float)),
                  **kw_data, label="Data")
    fill_plot(ax_e, bins_e, comp_bkg_E, "#FFA500")
    fill_plot(ax_e, bins_e, comp_Li9_E, "#00BFFF")
    fill_plot(ax_e, bins_e, comp_He8_E, "#9B59B6")
    step_plot(ax_e, bins_e, comp_bkg_E, **kw_bkg, label="Background")
    step_plot(ax_e, bins_e, comp_Li9_E, **kw_Li9, label=r"$^9$Li")
    step_plot(ax_e, bins_e, comp_He8_E, **kw_He8, label=r"$^8$He")
    step_plot(ax_e, bins_e, total_E,    **kw_sum, label="Total fit")
    ax_e.set_xlabel(r"$E_p$ (MeV)")
    ax_e.set_ylabel(r"Entries")
    ax_e.set_xlim(bins_e[0], bins_e[-1])
    ax_e.set_ylim(bottom=0)
    ax_e.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_e.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax_e.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
    ax_e.legend(loc="upper right")
    fig_e.tight_layout()

    # ── Delta-t plot ──────────────────────────────────────────────────────────
    fig_dt, ax_dt = plt.subplots(figsize=(7, 6))
    bins_dt_masked = np.r_[dt_centres - dt_width / 2,
                            dt_centres[-1] + dt_width[-1] / 2]
    ax_dt.errorbar(dt_centres, dt_sig, yerr=np.sqrt(dt_sig.astype(float)),
                   **kw_data, label="Data")
    fill_plot(ax_dt, bins_dt_masked, comp_bkg_dt, "#FFA500")
    fill_plot(ax_dt, bins_dt_masked, comp_Li9_dt, "#00BFFF")
    fill_plot(ax_dt, bins_dt_masked, comp_He8_dt, "#9B59B6")
    step_plot(ax_dt, bins_dt_masked, comp_bkg_dt, **kw_bkg, label="Background")
    step_plot(ax_dt, bins_dt_masked, comp_Li9_dt, **kw_Li9, label=r"$^9$Li")
    step_plot(ax_dt, bins_dt_masked, comp_He8_dt, **kw_He8, label=r"$^8$He")
    step_plot(ax_dt, bins_dt_masked, total_dt,    **kw_sum, label="Total fit")
    ax_dt.set_xlabel(r"$\Delta t_{\mu \to p}$ (s)")
    ax_dt.set_ylabel(r"Entries")
    ax_dt.set_xlim(bins_dt_masked[0], bins_dt_masked[-1])
    ax_dt.set_ylim(bottom=0)
    ax_dt.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_dt.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax_dt.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
    ax_dt.legend(loc="upper right")
    fig_dt.tight_layout()

    # ── 1D chi2 scans ─────────────────────────────────────────────────────────
    PARAM_NAMES = ["f_bkg", "f_Li9", "tau_Li9", "tau_He8"]
    LABELS = {
        "f_bkg":   r"$f_{\rm bkg}$",
        "f_Li9":   r"$f_{^9{\rm Li}}$",
        "tau_Li9": r"$\tau_{^9{\rm Li}}$ (ms)",
        "tau_He8": r"$\tau_{^8{\rm He}}$ (ms)",
    }
    SCALE = {"f_bkg": 1, "f_Li9": 1, "tau_Li9": 1e3, "tau_He8": 1e3}

    chi2_min  = chi2_func(result.x[np.isin(PARAM_NAMES, free_names)])
    best_free = result.x[[PARAM_NAMES.index(n) for n in free_names]]

    n_free  = len(free_names)
    n_cols  = min(n_free, 2)
    n_rows  = (n_free + 1) // 2
    fig_1d, axes_1d = plt.subplots(n_rows, n_cols,
                                   figsize=(6 * n_cols, 4 * n_rows),
                                   squeeze=False)

    for idx, name in enumerate(free_names):
        ax  = axes_1d[idx // n_cols][idx % n_cols]
        sc  = SCALE[name]
        p0  = result.x[PARAM_NAMES.index(name)]
        err = errors[PARAM_NAMES.index(name)]

        # Scan ± 4 sigma around best fit
        p_scan = np.linspace(p0 - 4 * err, p0 + 4 * err, 100)
        chi2_scan = np.empty_like(p_scan)

        for k, pval in enumerate(p_scan):
            # Build free_vals with this parameter replaced
            free_vals = best_free.copy()
            free_vals[idx] = pval
            chi2_scan[k] = chi2_func(free_vals) - chi2_min

        ax.plot(p_scan * sc, chi2_scan, color="black", lw=1.5)
        ax.axvline(p0 * sc,          color="black", lw=1.0, linestyle="--", alpha=0.6)
        ax.axvline((p0 - err) * sc,  color="#e05c5c", lw=1.0, linestyle=":")
        ax.axvline((p0 + err) * sc,  color="#e05c5c", lw=1.0, linestyle=":")
        ax.axhline(1.0, color="#e05c5c", lw=0.8, linestyle="--", alpha=0.5,
                   label=r"$\Delta\chi^2 = 1$")
        ax.axhline(4.0, color="#f0a500", lw=0.8, linestyle="--", alpha=0.5,
                   label=r"$\Delta\chi^2 = 4$")
        ax.set_xlabel(LABELS[name])
        ax.set_ylabel(r"$\Delta\chi^2$")
        ax.set_ylim(bottom=0)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
        ax.legend(fontsize=14, loc="upper center")

    # Hide unused panels if n_free is odd
    if n_free % 2 != 0 and n_cols == 2:
        axes_1d[-1][-1].set_visible(False)

    fig_1d.tight_layout()

    # ── 2D contour plots ──────────────────────────────────────────────────────
    pairs = [(i, j) for i in range(n_free) for j in range(n_free) if i < j]
    n_pairs = len(pairs)

    if n_pairs > 0:
        n_cols_2d = min(n_pairs, 3)
        n_rows_2d = (n_pairs + n_cols_2d - 1) // n_cols_2d
        fig_2d, axes_2d = plt.subplots(n_rows_2d, n_cols_2d,
                                       figsize=(5 * n_cols_2d, 4.5 * n_rows_2d),
                                       squeeze=False)

        for plot_idx, (i, j) in enumerate(pairs):
            ax   = axes_2d[plot_idx // n_cols_2d][plot_idx % n_cols_2d]
            ni, nj = free_names[i], free_names[j]
            p0i  = result.x[PARAM_NAMES.index(ni)]
            p0j  = result.x[PARAM_NAMES.index(nj)]
            erri = errors[PARAM_NAMES.index(ni)]
            errj = errors[PARAM_NAMES.index(nj)]
            sci  = SCALE[ni]
            scj  = SCALE[nj]

            # Grid: ± 3 sigma on each axis, 40x40
            pi_vals = np.linspace(p0i - 3 * erri, p0i + 3 * erri, 40)
            pj_vals = np.linspace(p0j - 3 * errj, p0j + 3 * errj, 40)
            PI, PJ  = np.meshgrid(pi_vals, pj_vals)
            Z       = np.empty_like(PI)

            for ki in range(PI.shape[0]):
                for kj in range(PI.shape[1]):
                    free_vals      = best_free.copy()
                    free_vals[i]   = PI[ki, kj]
                    free_vals[j]   = PJ[ki, kj]
                    Z[ki, kj]      = chi2_func(free_vals) - chi2_min

            # 1-sigma and 2-sigma contours for 2 DOF: Delta chi2 = 2.30, 6.18
            contour_f = ax.contourf(PI * sci, PJ * scj, Z,
                                    levels=[0, 2.30, 6.18],
                                    colors=["#4a90d9", "#a8c8f0"],
                                    alpha=0.4)
            contour_l = ax.contour(PI * sci, PJ * scj, Z,
                                   levels=[2.30, 6.18],
                                   colors=["#1a5fa8", "#4a90d9"],
                                   linewidths=1.5)
            ax.clabel(contour_l,
                      fmt={2.30: r"$1\sigma$", 6.18: r"$2\sigma$"},
                      fontsize=13)
            ax.plot(p0i * sci, p0j * scj, "k+", ms=10, mew=1.5,
                    label="Best fit")

            ax.set_xlabel(LABELS[ni])
            ax.set_ylabel(LABELS[nj])
            ax.xaxis.set_minor_locator(AutoMinorLocator(5))
            ax.yaxis.set_minor_locator(AutoMinorLocator(5))
            ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
            ax.legend(fontsize=14, loc="upper right")

        # Hide unused panels
        for plot_idx in range(n_pairs, n_rows_2d * n_cols_2d):
            axes_2d[plot_idx // n_cols_2d][plot_idx % n_cols_2d].set_visible(False)

        fig_2d.tight_layout()

    plt.show()

if __name__ == "__main__":
    set_latex_style()

    args = parse_args()

    fin = uproot.open(args.data)

    tin = fin["background_events"]
    din = tin.arrays(library="np")
    
    background_run_id = din["run_id"]
    background_posx_p = din["posx_p"]
    background_posy_p = din["posy_p"]
    background_posz_p = din["posz_p"]
    background_sec_p = din["sec_p"]
    background_nsec_p = din["nsec_p"]
    background_e_p = din["e_p"]
    background_dlat_mu2p = din["dlat_mu2p"]
    background_dt_mu2p = din["dt_mu2p"]
    background_posx_d = din["posx_d"]
    background_posy_d = din["posy_d"]
    background_posz_d = din["posz_d"]
    background_sec_d = din["sec_d"]
    background_nsec_d = din["nsec_d"]
    background_e_d = din["e_d"]
    background_dlat_mu2d = din["dlat_mu2d"]
    background_dt_mu2d = din["dt_mu2d"]

    tin = fin["signal_events"]
    din = tin.arrays(library="np")

    signal_run_id = din["run_id"]
    signal_posx_p = din["posx_p"]
    signal_posy_p = din["posy_p"]
    signal_posz_p = din["posz_p"]
    signal_sec_p = din["sec_p"]
    signal_nsec_p = din["nsec_p"]
    signal_e_p = din["e_p"]
    signal_dlat_mu2p = din["dlat_mu2p"]
    signal_dt_mu2p = din["dt_mu2p"]
    signal_posx_d = din["posx_d"]
    signal_posy_d = din["posy_d"]
    signal_posz_d = din["posz_d"]
    signal_sec_d = din["sec_d"]
    signal_nsec_d = din["nsec_d"]
    signal_e_d = din["e_d"]
    signal_dlat_mu2d = din["dlat_mu2d"]
    signal_dt_mu2d = din["dt_mu2d"]

    fin = uproot.open(args.simulation)

    tin = fin["cosmogenics"]
    din = tin.arrays(library="np")

    mc_posx_p = din["posx_p"]
    mc_posy_p = din["posy_p"]
    mc_posz_p = din["posz_p"]
    mc_e_p = din["e_p"]
    mc_element = din["element"]

    # bins = nmo_analysis_bins()
    bins = np.linspace(0.0, 12.0, 51)
    histmax = 0.0

    fig, ax = plt.subplots(figsize=(7, 6))

    signal_hist, edges = np.histogram(signal_e_p, bins=bins)
    histmax = np.max([histmax, np.max(signal_hist)])
    ax.fill_between(bins, np.r_[signal_hist, signal_hist[-1]], step="post", color="#648fff", alpha=0.075, zorder=1)
    ax.step(bins, np.r_[signal_hist, signal_hist[-1]], where="post", color="#648fff", linestyle="--", linewidth=1.5, zorder=2, label="Cosmogenic enriched region")
    background_hist, edges = np.histogram(background_e_p, bins=bins)
    histmax = np.max([histmax, np.max(background_hist)])
    ax.fill_between(bins, np.r_[background_hist, background_hist[-1]], step="post", color="#ff6464", alpha=0.075, zorder=1)
    ax.step(bins, np.r_[background_hist, background_hist[-1]], where="post", color="#ff6464", linestyle="--", linewidth=1.5, zorder=2, label="Cosmogenic depleted region")
    hist_difference = signal_hist - background_hist
    ax.step(bins, np.r_[hist_difference, hist_difference[-1]], where="post", color="#000000", linestyle="-", linewidth=1.5, zorder=2, label="Difference")

    ax.set_xlabel(r"$E_{p}$ (MeV)")
    ax.set_ylabel(r"Entries")
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(bottom=0.0, top=histmax * 1.25)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.legend(loc="upper right")

    fig.tight_layout()
    fig.show()

    fig, ax = plt.subplots(figsize=(7, 6))

    li9_mask = mc_element == "Li9"
    he8_mask = mc_element == "He8"

    li9_hist, edeges = np.histogram(mc_e_p[li9_mask], bins=bins)
    histmax = np.max([histmax, np.max(li9_hist)])
    ax.fill_between(bins, np.r_[li9_hist, li9_hist[-1]], step="post", color="#ff6464", alpha=0.075, zorder=1)
    ax.step(bins, np.r_[li9_hist, li9_hist[-1]], where="post", color="#ff6464", linestyle="-", linewidth=1.5, zorder=2, label=r"$^9$Li")
    he8_hist, edges = np.histogram(mc_e_p[he8_mask], bins=bins)
    histmax = np.max([histmax, np.max(he8_hist)])
    ax.fill_between(bins, np.r_[he8_hist, he8_hist[-1]], step="post", color="#84d040", alpha=0.075, zorder=1)
    ax.step(bins, np.r_[he8_hist, he8_hist[-1]], where="post", color="#84d040", linestyle="-", linewidth=1.5, zorder=2, label=r"$^8$He")

    ax.set_xlabel(r"$E_{p}$ (MeV)")
    ax.set_ylabel(r"Entries")
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(bottom=0.0, top=histmax * 1.25)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.legend(loc="upper right")

    fig.tight_layout()
    fig.show()

    bins = np.linspace(0.0, 1.2, 101)
    histmax = 0.0

    fig, ax = plt.subplots(figsize=(7, 6))

    signal_hist = np.histogram(signal_dt_mu2p, bins=bins)[0]
    histmax = np.max([histmax, np.max(signal_hist)])
    ax.fill_between(bins, np.r_[signal_hist, signal_hist[-1]], step="post", color="#648fff", alpha=0.075, zorder=1)
    ax.step(bins, np.r_[signal_hist, signal_hist[-1]], where="post", color="#648fff", linestyle="--", linewidth=1.5, zorder=2, label="Cosmogenic enriched region")
    background_hist = np.histogram(background_dt_mu2p, bins=bins)[0]
    histmax = np.max([histmax, np.max(background_hist)])
    ax.fill_between(bins, np.r_[background_hist, background_hist[-1]], step="post", color="#ff6464", alpha=0.075, zorder=1)
    ax.step(bins, np.r_[background_hist, background_hist[-1]], where="post", color="#ff6464", linestyle="--", linewidth=1.5, zorder=2, label="Cosmogenic depleted region")

    ax.set_xlabel(r"$\Delta t$ (s)")
    ax.set_ylabel(r"Entries")
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(bottom=0.0, top=histmax * 1.25)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.legend(loc="upper right")

    fig.tight_layout()
    fig.show()

    # S = f_{bkg} * S_{bkg} + (1 - f_{bkg}) [f_{Li9} S_{Li9} + (1 - f_{Li9}) * S_{He8}]
    # dt = f_{bkg} * dt_{bkg} + (1 - f_{bkg}) [f_{Li9} * exp(-t/\tau_{Li9}) + (1 - f_{Li9}) * exp(-t/\tau_{He8})]

    bins_e  = np.linspace(0.0, 12.0, 51)
    bins_dt = np.linspace(0.0, 2.0, 101)

    result, errors, S_sig, S_bkg, S_Li9, S_He8, dt_sig, chi2_func, free_names = run_fit(
        signal_e_p, signal_dt_mu2p,
        background_e_p, mc_e_p, mc_element,
        bins_e, 
        background_dt_mu2p, bins_dt,
        # fixed={"tau_Li9": 0.257, "tau_He8": 0.172}
    )

    plot_fit(
        result, errors, S_sig, S_bkg, S_Li9, S_He8, dt_sig,
        bins_e, bins_dt, chi2_func, free_names,
        # fixed={"tau_Li9": 0.257, "tau_He8": 0.172}
    )

    plt.show()