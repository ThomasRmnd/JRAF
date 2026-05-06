"""
cosmogenic_fit.py
=================
Chi2 fit of cosmogenic isotopes (Li9 / He8) in a reactor antineutrino detector.

Model
-----
    S(E) = f_bkg * S_bkg(E)
         + (1 - f_bkg) * [f_Li9 * S_Li9(E) + (1 - f_Li9) * S_He8(E)]

    dt(t) = f_bkg * Uniform(t)
          + (1 - f_bkg) * [f_Li9 * Exp(t; tau_Li9) + (1 - f_Li9) * Exp(t; tau_He8)]

Free parameters
---------------
    f_bkg   : fraction of background in signal region         in [0, 1]
    f_Li9   : fraction of Li9 among cosmogenic isotopes       in [0, 1]
    tau_Li9 : mean lifetime of Li9                            in seconds
    tau_He8 : mean lifetime of He8                            in seconds

Any subset of these can be fixed at the call site via fixed={...}.

Usage
-----
    python cosmogenic_fit.py --data data.root --simulation simu.root
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numdifftools as nd
import numpy as np
from scipy.linalg import inv
from scipy.optimize import minimize, OptimizeResult
import uproot


# ══════════════════════════════════════════════════════════════════════════════
# Constants & PDG reference values
# ══════════════════════════════════════════════════════════════════════════════

# Canonical parameter order used throughout
PARAM_NAMES: tuple[str, ...] = ("f_bkg", "f_Li9", "tau_Li9", "tau_He8")

PDG_LIFETIMES: dict[str, float] = {
    "tau_Li9": 0.257,   # s  (T_{1/2} = 178 ms / ln2)
    "tau_He8": 0.172,   # s  (T_{1/2} = 119 ms / ln2)
}

DEFAULT_INITIAL_VALUES: dict[str, float] = {
    "f_bkg":   0.08,
    "f_Li9":   0.90,
    "tau_Li9": PDG_LIFETIMES["tau_Li9"],
    "tau_He8": PDG_LIFETIMES["tau_He8"],
}

PARAM_BOUNDS: dict[str, tuple[float, float]] = {
    "f_bkg":   (0.0,  1.0),
    "f_Li9":   (0.0,  1.0),
    "tau_Li9": (0.05, 1.0),
    "tau_He8": (0.05, 1.0),
}

PARAM_LABELS: dict[str, str] = {
    "f_bkg":   r"$f_{\rm bkg}$",
    "f_Li9":   r"$f_{^9{\rm Li}}$",
    "tau_Li9": r"$\tau_{^9{\rm Li}}$ (ms)",
    "tau_He8": r"$\tau_{^8{\rm He}}$ (ms)",
}

# Display scale: internal unit → plot unit
PARAM_SCALE: dict[str, float] = {
    "f_bkg":   1.0,
    "f_Li9":   1.0,
    "tau_Li9": 1e3,   # s → ms
    "tau_He8": 1e3,
}

PARAM_FMT: dict[str, tuple[str, str]] = {
    "f_bkg":   (".4f", ""),
    "f_Li9":   (".4f", ""),
    "tau_Li9": (".1f", " ms"),
    "tau_He8": (".1f", " ms"),
}

# Minimum dt used in the fit (remove prompt coincidences)
DT_FIT_MIN: float = 0.060   # seconds

# Delta-chi2 thresholds for confidence regions
CONTOUR_1SIGMA_2DOF: float = 2.30
CONTOUR_2SIGMA_2DOF: float = 6.18


# ══════════════════════════════════════════════════════════════════════════════
# Matplotlib style
# ══════════════════════════════════════════════════════════════════════════════

def set_latex_style() -> None:
    """Apply publication-quality LaTeX matplotlib style."""
    mpl.rcParams.update({
        "text.usetex":        True,
        "font.family":        "serif",
        "font.serif":         ["Computer Modern Serif"],
        "mathtext.fontset":   "cm",
        "font.size":          22,
        "axes.labelsize":     22,
        "axes.titlesize":     22,
        "xtick.labelsize":    18,
        "ytick.labelsize":    18,
        "legend.fontsize":    18,
        "axes.linewidth":     1.35,
        "xtick.direction":    "in",
        "ytick.direction":    "in",
        "xtick.major.size":   10,
        "ytick.major.size":   10,
        "xtick.minor.size":   5,
        "ytick.minor.size":   5,
        "xtick.major.width":  1.25,
        "ytick.major.width":  1.25,
        "xtick.minor.width":  0.75,
        "ytick.minor.width":  0.75,
        "xtick.top":          True,
        "ytick.right":        True,
        "legend.frameon":     False,
        "figure.figsize":     (8, 6),
        "figure.dpi":         120,
        "savefig.bbox":       "tight",
        "savefig.dpi":        300,
    })


# ══════════════════════════════════════════════════════════════════════════════
# Data containers
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class CosmogenicData:
    """Raw arrays loaded from the data ROOT file."""
    # Prompt quantities
    run_id:      np.ndarray
    posx_p:      np.ndarray
    posy_p:      np.ndarray
    posz_p:      np.ndarray
    sec_p:       np.ndarray
    nsec_p:      np.ndarray
    e_p:         np.ndarray
    dlat_mu2p:   np.ndarray
    dt_mu2p:     np.ndarray
    # Delayed quantities
    posx_d:      np.ndarray
    posy_d:      np.ndarray
    posz_d:      np.ndarray
    sec_d:       np.ndarray
    nsec_d:      np.ndarray
    e_d:         np.ndarray
    dlat_mu2d:   np.ndarray
    dt_mu2d:     np.ndarray

    @classmethod
    def from_tree(cls, tree) -> "CosmogenicData":
        d = tree.arrays(library="np")
        return cls(**{f: d[f] for f in cls.__dataclass_fields__})


@dataclass
class SimulationData:
    """Raw arrays loaded from the simulation ROOT file."""
    posx_p:   np.ndarray
    posy_p:   np.ndarray
    posz_p:   np.ndarray
    e_p:      np.ndarray
    element:  np.ndarray

    @classmethod
    def from_tree(cls, tree) -> "SimulationData":
        d = tree.arrays(library="np")
        return cls(**{f: d[f] for f in cls.__dataclass_fields__})


@dataclass
class FitTemplates:
    """
    Pre-histogrammed templates and masked dt arrays used by the chi2.
    Built once; reused across minimisation calls.
    """
    # Energy templates (raw counts)
    S_sig:   np.ndarray
    S_bkg:   np.ndarray
    S_Li9:   np.ndarray
    S_He8:   np.ndarray
    # Normalised energy templates
    S_bkg_n: np.ndarray
    S_Li9_n: np.ndarray
    S_He8_n: np.ndarray
    # Total signal counts (fixed during fit)
    N_sig:   int
    # dt arrays after DT_FIT_MIN cut
    dt_sig:      np.ndarray
    dt_centres:  np.ndarray
    dt_width:    np.ndarray

    @staticmethod
    def _safe_norm(h: np.ndarray) -> np.ndarray:
        n = h.sum()
        return h / n if n > 0 else h.astype(float)

    @classmethod
    def build(
        cls,
        signal:     CosmogenicData,
        background: CosmogenicData,
        mc:         SimulationData,
        bins_e:     np.ndarray,
        bins_dt:    np.ndarray,
    ) -> "FitTemplates":
        S_sig, _ = np.histogram(signal.e_p,     bins=bins_e)
        S_bkg, _ = np.histogram(background.e_p, bins=bins_e)
        S_Li9, _ = np.histogram(mc.e_p[mc.element == "Li9"], bins=bins_e)
        S_He8, _ = np.histogram(mc.e_p[mc.element == "He8"], bins=bins_e)

        dt_centres_full = 0.5 * (bins_dt[:-1] + bins_dt[1:])
        dt_width_full   = np.diff(bins_dt)
        dt_mask         = dt_centres_full > DT_FIT_MIN

        dt_sig_full, _ = np.histogram(signal.dt_mu2p, bins=bins_dt)

        return cls(
            S_sig   = S_sig,
            S_bkg   = S_bkg,
            S_Li9   = S_Li9,
            S_He8   = S_He8,
            S_bkg_n = cls._safe_norm(S_bkg.astype(float)),
            S_Li9_n = cls._safe_norm(S_Li9.astype(float)),
            S_He8_n = cls._safe_norm(S_He8.astype(float)),
            N_sig   = int(S_sig.sum()),
            dt_sig      = dt_sig_full[dt_mask],
            dt_centres  = dt_centres_full[dt_mask],
            dt_width    = dt_width_full[dt_mask],
        )


@dataclass
class FitResult:
    """Full fit output: best-fit values, uncertainties, and diagnostics."""
    # Best-fit values in canonical PARAM_NAMES order
    params:     dict[str, float]
    errors:     dict[str, float]
    fixed:      dict[str, float]
    free_names: list[str]
    chi2_min:   float
    ndof:       int
    converged:  bool
    # The callable chi2(free_vals) for post-fit scans / contours
    chi2_func:  object = field(repr=False)

    @property
    def chi2_per_ndof(self) -> float:
        return self.chi2_min / self.ndof if self.ndof > 0 else float("nan")

    def summary(self) -> str:
        lines = [
            f"Fit converged : {self.converged}",
            f"chi2 / ndof   = {self.chi2_min:.1f} / {self.ndof}"
            f" = {self.chi2_per_ndof:.3f}",
        ]
        for name in PARAM_NAMES:
            fmt, unit = PARAM_FMT[name]
            sc        = PARAM_SCALE[name]
            val       = self.params[name] * sc
            if name in self.fixed:
                lines.append(f"  {name:8s} = {val:{fmt}}{unit}  [fixed]")
            else:
                err = self.errors[name] * sc
                pdg = f"  (PDG: {PDG_LIFETIMES[name]*1e3:.0f} ms)" \
                      if name in PDG_LIFETIMES else ""
                lines.append(
                    f"  {name:8s} = {val:{fmt}}{unit}"
                    f"  +/-  {err:{fmt}}{unit}{pdg}"
                )
        return "\n".join(lines)


# ══════════════════════════════════════════════════════════════════════════════
# Chi2 construction
# ══════════════════════════════════════════════════════════════════════════════

class Chi2Builder:
    """
    Builds the chi2 callable from pre-computed templates.

    The chi2 is a function of the *free* parameters only (a 1-D numpy array
    in the order given by ``free_names``).  Fixed parameters are baked in via
    closure so the minimiser never sees them.
    """

    def __init__(self, templates: FitTemplates, fixed: dict[str, float]) -> None:
        self.templates  = templates
        self.fixed      = fixed
        self.free_names = [p for p in PARAM_NAMES if p not in fixed]

    @property
    def free_bounds(self) -> list[tuple[float, float]]:
        return [PARAM_BOUNDS[p] for p in self.free_names]

    @property
    def free_x0(self) -> list[float]:
        return [DEFAULT_INITIAL_VALUES[p] for p in self.free_names]

    def _unpack(self, free_vals: np.ndarray) -> tuple[float, float, float, float]:
        d = dict(self.fixed)
        d.update(dict(zip(self.free_names, free_vals)))
        return d["f_bkg"], d["f_Li9"], d["tau_Li9"], d["tau_He8"]

    @staticmethod
    def _safe_norm_arr(a: np.ndarray) -> np.ndarray:
        s = a.sum()
        return a / s if s > 0 else a

    def __call__(self, free_vals: np.ndarray) -> float:
        t = self.templates
        f_bkg, f_Li9, tau_Li9, tau_He8 = self._unpack(free_vals)

        # ── Energy chi2 ───────────────────────────────────────────────────────
        E_model = t.N_sig * (
            f_bkg       * t.S_bkg_n
            + (1-f_bkg) * (f_Li9 * t.S_Li9_n + (1-f_Li9) * t.S_He8_n)
        )
        sigma_E = np.sqrt(
            f_bkg       * t.S_bkg.astype(float)
            + (1-f_bkg) * (f_Li9 * t.S_Li9.astype(float)
                           + (1-f_Li9) * t.S_He8.astype(float))
        )
        sigma_E = np.where(sigma_E > 0, sigma_E, 1.0)
        chi2_E  = np.sum(((t.S_sig - E_model) / sigma_E) ** 2)

        # ── dt chi2 ───────────────────────────────────────────────────────────
        sn = self._safe_norm_arr

        exp_Li9 = np.exp(-t.dt_centres / tau_Li9) * t.dt_width
        exp_He8 = np.exp(-t.dt_centres / tau_He8) * t.dt_width
        flat    = np.ones_like(t.dt_centres)       * t.dt_width

        dt_model = t.N_sig * (
            f_bkg       * sn(flat)
            + (1-f_bkg) * (f_Li9 * sn(exp_Li9) + (1-f_Li9) * sn(exp_He8))
        )
        sigma_dt_shape = t.N_sig * (
            f_bkg       * sn(flat)
            + (1-f_bkg) * (f_Li9 * sn(exp_Li9) + (1-f_Li9) * sn(exp_He8))
        )
        sigma_dt = np.where(sigma_dt_shape > 0, np.sqrt(sigma_dt_shape), 1.0)
        chi2_dt  = np.sum(((t.dt_sig - dt_model) / sigma_dt) ** 2)

        return float(chi2_E + chi2_dt)


# ══════════════════════════════════════════════════════════════════════════════
# Fitter
# ══════════════════════════════════════════════════════════════════════════════

class CosmogenicFitter:
    """
    Wraps Chi2Builder + scipy.minimize + Hessian uncertainty estimation.

    Parameters
    ----------
    templates : FitTemplates
    fixed     : dict mapping parameter name → fixed value (optional)
    """

    def __init__(
        self,
        templates: FitTemplates,
        fixed: Optional[dict[str, float]] = None,
    ) -> None:
        self.templates = templates
        self.fixed     = fixed or {}
        self.builder   = Chi2Builder(templates, self.fixed)

    def fit(self) -> FitResult:
        chi2 = self.builder

        opt: OptimizeResult = minimize(
            chi2,
            chi2.free_x0,
            method  = "L-BFGS-B",
            bounds  = chi2.free_bounds,
            options = {"ftol": 1e-12, "gtol": 1e-8, "maxiter": 100_000},
        )

        # ── Uncertainties from Hessian of chi2 ────────────────────────────────
        H      = nd.Hessian(chi2)(opt.x)
        cov    = 2.0 * inv(H)
        errors = np.sqrt(np.diag(cov))

        # ── Pack into canonical dicts ─────────────────────────────────────────
        free_params = dict(zip(chi2.free_names, opt.x))
        free_errors = dict(zip(chi2.free_names, errors))

        params = {**self.fixed, **free_params}
        errs   = {**{n: 0.0 for n in self.fixed}, **free_errors}

        ndof = (
            len(self.templates.S_sig)
            + len(self.templates.dt_sig)
            - len(chi2.free_names)
        )

        return FitResult(
            params     = params,
            errors     = errs,
            fixed      = self.fixed,
            free_names = chi2.free_names,
            chi2_min   = float(opt.fun),
            ndof       = ndof,
            converged  = bool(opt.success),
            chi2_func  = chi2,
        )


# ══════════════════════════════════════════════════════════════════════════════
# Plotting helpers
# ══════════════════════════════════════════════════════════════════════════════

def _step_plot(ax, bins: np.ndarray, values: np.ndarray, **kwargs) -> None:
    ax.step(bins, np.r_[values, values[-1]], where="post", **kwargs)


def _fill_plot(ax, bins: np.ndarray, values: np.ndarray, color: str) -> None:
    ax.fill_between(bins, np.r_[values, values[-1]],
                    step="post", color=color, alpha=0.15)


def _decorate_ax(ax, xlabel: str, ylabel: str) -> None:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.legend(loc="upper right")


# ══════════════════════════════════════════════════════════════════════════════
# Fit result plotter
# ══════════════════════════════════════════════════════════════════════════════

class FitPlotter:
    """
    Produces all diagnostic plots from a FitResult + FitTemplates.

    Plots produced
    --------------
    1. Energy spectrum with fitted components
    2. dt spectrum with fitted components
    3. 1-D chi2 scans for each free parameter
    4. 2-D chi2 contour plots for each pair of free parameters
    """

    # Colour palette
    COLOR_BKG  = "#FFA500"
    COLOR_LI9  = "#00BFFF"
    COLOR_HE8  = "#9B59B6"
    COLOR_SUM  = "black"
    COLOR_DATA = "black"

    def __init__(
        self,
        result:    FitResult,
        templates: FitTemplates,
        bins_e:    np.ndarray,
        bins_dt:   np.ndarray,
    ) -> None:
        self.result    = result
        self.templates = templates
        self.bins_e    = bins_e
        self.bins_dt   = bins_dt
        self._build_model_components()

    # ── Internal model reconstruction ─────────────────────────────────────────

    @staticmethod
    def _safe_norm(h: np.ndarray) -> np.ndarray:
        n = h.sum()
        return h / n if n > 0 else h.astype(float)

    @staticmethod
    def _safe_norm_arr(a: np.ndarray) -> np.ndarray:
        s = a.sum()
        return a / s if s > 0 else a

    def _build_model_components(self) -> None:
        p  = self.result.params
        t  = self.templates
        sn = self._safe_norm

        f_bkg, f_Li9 = p["f_bkg"], p["f_Li9"]
        tau_Li9, tau_He8 = p["tau_Li9"], p["tau_He8"]
        N = t.N_sig

        # Energy components
        self.comp_bkg_E = N * f_bkg           * sn(t.S_bkg.astype(float))
        self.comp_Li9_E = N * (1-f_bkg)*f_Li9 * sn(t.S_Li9.astype(float))
        self.comp_He8_E = N * (1-f_bkg)*(1-f_Li9) * sn(t.S_He8.astype(float))
        self.total_E    = self.comp_bkg_E + self.comp_Li9_E + self.comp_He8_E

        # dt components (masked)
        sna = self._safe_norm_arr
        dc, dw = t.dt_centres, t.dt_width

        exp_Li9 = np.exp(-dc / tau_Li9) * dw
        exp_He8 = np.exp(-dc / tau_He8) * dw
        flat    = np.ones_like(dc) * dw

        self.comp_bkg_dt = N * f_bkg           * sna(flat)
        self.comp_Li9_dt = N * (1-f_bkg)*f_Li9 * sna(exp_Li9)
        self.comp_He8_dt = N * (1-f_bkg)*(1-f_Li9) * sna(exp_He8)
        self.total_dt    = (self.comp_bkg_dt + self.comp_Li9_dt
                            + self.comp_He8_dt)

        # Masked dt bins for step plots
        self.bins_dt_masked = np.r_[
            dc - dw/2,
            dc[-1] + dw[-1]/2
        ]

    # ── Public plot methods ────────────────────────────────────────────────────

    def plot_energy(self) -> plt.Figure:
        t = self.templates
        E_centres = 0.5 * (self.bins_e[:-1] + self.bins_e[1:])

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.errorbar(E_centres, t.S_sig,
                    yerr=np.sqrt(t.S_sig.astype(float)),
                    fmt="ko", ms=4, capsize=3, lw=1.0, zorder=5, label="Data")

        _fill_plot(ax, self.bins_e, self.comp_bkg_E, self.COLOR_BKG)
        _fill_plot(ax, self.bins_e, self.comp_Li9_E, self.COLOR_LI9)
        _fill_plot(ax, self.bins_e, self.comp_He8_E, self.COLOR_HE8)
        _step_plot(ax, self.bins_e, self.comp_bkg_E,
                   color=self.COLOR_BKG, lw=1.5, label="Background")
        _step_plot(ax, self.bins_e, self.comp_Li9_E,
                   color=self.COLOR_LI9, lw=1.5, label=r"$^9$Li")
        _step_plot(ax, self.bins_e, self.comp_He8_E,
                   color=self.COLOR_HE8, lw=1.5, label=r"$^8$He")
        _step_plot(ax, self.bins_e, self.total_E,
                   color=self.COLOR_SUM, lw=1.5, linestyle="dotted",
                   label="Total fit")

        ax.set_xlim(self.bins_e[0], self.bins_e[-1])
        ax.set_ylim(bottom=0)
        _decorate_ax(ax, r"$E_p$ (MeV)", r"Entries")
        fig.tight_layout()
        return fig

    def plot_dt(self) -> plt.Figure:
        t  = self.templates
        dc = t.dt_centres

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.errorbar(dc, t.dt_sig,
                    yerr=np.sqrt(t.dt_sig.astype(float)),
                    fmt="ko", ms=4, capsize=3, lw=1.0, zorder=5, label="Data")

        bm = self.bins_dt_masked
        _fill_plot(ax, bm, self.comp_bkg_dt, self.COLOR_BKG)
        _fill_plot(ax, bm, self.comp_Li9_dt, self.COLOR_LI9)
        _fill_plot(ax, bm, self.comp_He8_dt, self.COLOR_HE8)
        _step_plot(ax, bm, self.comp_bkg_dt,
                   color=self.COLOR_BKG, lw=1.5, label="Background")
        _step_plot(ax, bm, self.comp_Li9_dt,
                   color=self.COLOR_LI9, lw=1.5, label=r"$^9$Li")
        _step_plot(ax, bm, self.comp_He8_dt,
                   color=self.COLOR_HE8, lw=1.5, label=r"$^8$He")
        _step_plot(ax, bm, self.total_dt,
                   color=self.COLOR_SUM, lw=1.5, linestyle="dotted",
                   label="Total fit")

        ax.set_xlim(bm[0], bm[-1])
        ax.set_ylim(bottom=0)
        _decorate_ax(ax, r"$\Delta t_{\mu \to p}$ (s)", r"Entries")
        fig.tight_layout()
        return fig

    def plot_chi2_scans(self) -> plt.Figure:
        """1-D Delta-chi2 profiles for each free parameter."""
        r          = self.result
        free_names = r.free_names
        n_free     = len(free_names)

        if n_free == 0:
            return None

        chi2_min  = r.chi2_min
        best_free = np.array([r.params[n] for n in free_names])

        n_cols = min(n_free, 2)
        n_rows = (n_free + 1) // 2
        fig, axes = plt.subplots(n_rows, n_cols,
                                 figsize=(6*n_cols, 4*n_rows),
                                 squeeze=False)

        for idx, name in enumerate(free_names):
            ax  = axes[idx // n_cols][idx % n_cols]
            sc  = PARAM_SCALE[name]
            p0  = r.params[name]
            err = r.errors[name]

            p_scan    = np.linspace(p0 - 4*err, p0 + 4*err, 120)
            chi2_scan = np.empty_like(p_scan)

            for k, pval in enumerate(p_scan):
                fv      = best_free.copy()
                fv[idx] = pval
                chi2_scan[k] = r.chi2_func(fv) - chi2_min

            ax.plot(p_scan * sc, chi2_scan, color="black", lw=1.5)
            ax.axvline(p0 * sc,          color="black",   lw=1.0,
                       linestyle="--", alpha=0.6)
            ax.axvline((p0 - err) * sc,  color="#e05c5c", lw=1.0,
                       linestyle=":")
            ax.axvline((p0 + err) * sc,  color="#e05c5c", lw=1.0,
                       linestyle=":")
            ax.axhline(1.0, color="#e05c5c", lw=0.8, linestyle="--",
                       alpha=0.5, label=r"$\Delta\chi^2 = 1$")
            ax.axhline(4.0, color="#f0a500", lw=0.8, linestyle="--",
                       alpha=0.5, label=r"$\Delta\chi^2 = 4$")
            ax.set_xlabel(PARAM_LABELS[name])
            ax.set_ylabel(r"$\Delta\chi^2$")
            ax.set_ylim(bottom=0)
            ax.xaxis.set_minor_locator(AutoMinorLocator(5))
            ax.yaxis.set_minor_locator(AutoMinorLocator(5))
            ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
            ax.legend(fontsize=14, loc="upper center")

        if n_free % 2 != 0 and n_cols == 2:
            axes[-1][-1].set_visible(False)

        fig.tight_layout()
        return fig

    def plot_contours(self) -> Optional[plt.Figure]:
        """2-D Delta-chi2 contour plots for each pair of free parameters."""
        r          = self.result
        free_names = r.free_names
        n_free     = len(free_names)
        pairs      = [(i, j) for i in range(n_free)
                              for j in range(n_free) if i < j]

        if not pairs:
            return None

        chi2_min  = r.chi2_min
        best_free = np.array([r.params[n] for n in free_names])

        n_pairs   = len(pairs)
        n_cols_2d = min(n_pairs, 3)
        n_rows_2d = (n_pairs + n_cols_2d - 1) // n_cols_2d

        fig, axes = plt.subplots(n_rows_2d, n_cols_2d,
                                 figsize=(5*n_cols_2d, 4.5*n_rows_2d),
                                 squeeze=False)

        for plot_idx, (i, j) in enumerate(pairs):
            ax = axes[plot_idx // n_cols_2d][plot_idx % n_cols_2d]
            ni, nj = free_names[i], free_names[j]

            p0i, p0j   = r.params[ni], r.params[nj]
            erri, errj = r.errors[ni], r.errors[nj]
            sci,  scj  = PARAM_SCALE[ni], PARAM_SCALE[nj]

            pi_vals = np.linspace(p0i - 3*erri, p0i + 3*erri, 40)
            pj_vals = np.linspace(p0j - 3*errj, p0j + 3*errj, 40)
            PI, PJ  = np.meshgrid(pi_vals, pj_vals)
            Z       = np.empty_like(PI)

            for ki in range(PI.shape[0]):
                for kj in range(PI.shape[1]):
                    fv    = best_free.copy()
                    fv[i] = PI[ki, kj]
                    fv[j] = PJ[ki, kj]
                    Z[ki, kj] = r.chi2_func(fv) - chi2_min

            cf = ax.contourf(PI*sci, PJ*scj, Z,
                             levels=[0, CONTOUR_1SIGMA_2DOF, CONTOUR_2SIGMA_2DOF],
                             colors=["#4a90d9", "#a8c8f0"], alpha=0.4)
            cl = ax.contour(PI*sci, PJ*scj, Z,
                            levels=[CONTOUR_1SIGMA_2DOF, CONTOUR_2SIGMA_2DOF],
                            colors=["#1a5fa8", "#4a90d9"], linewidths=1.5)
            ax.clabel(cl, fmt={CONTOUR_1SIGMA_2DOF: r"$1\sigma$",
                               CONTOUR_2SIGMA_2DOF: r"$2\sigma$"},
                      fontsize=13)
            ax.plot(p0i*sci, p0j*scj, "k+", ms=10, mew=1.5,
                    label="Best fit")
            ax.set_xlabel(PARAM_LABELS[ni])
            ax.set_ylabel(PARAM_LABELS[nj])
            ax.xaxis.set_minor_locator(AutoMinorLocator(5))
            ax.yaxis.set_minor_locator(AutoMinorLocator(5))
            ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
            ax.legend(fontsize=14, loc="upper right")

        for plot_idx in range(n_pairs, n_rows_2d * n_cols_2d):
            axes[plot_idx // n_cols_2d][plot_idx % n_cols_2d].set_visible(False)

        fig.tight_layout()
        return fig

    def plot_all(self) -> None:
        """Render all four diagnostic figures."""
        self.plot_energy()
        self.plot_dt()
        self.plot_chi2_scans()
        self.plot_contours()
        plt.show()


# ══════════════════════════════════════════════════════════════════════════════
# Exploratory pre-fit plots
# ══════════════════════════════════════════════════════════════════════════════

def plot_raw_spectra(
    signal:     CosmogenicData,
    background: CosmogenicData,
    mc:         SimulationData,
    bins_e:     np.ndarray,
    bins_dt:    np.ndarray,
) -> None:
    """Three overview plots produced before the fit."""

    # ── 1. Energy: signal vs background ──────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 6))
    histmax = 0.0

    for data, color, label in [
        (signal.e_p,     "#648fff", "Cosmogenic enriched region"),
        (background.e_p, "#ff6464", "Cosmogenic depleted region"),
    ]:
        h, _ = np.histogram(data, bins=bins_e)
        histmax = max(histmax, h.max())
        ax.fill_between(bins_e, np.r_[h, h[-1]], step="post",
                        color=color, alpha=0.075, zorder=1)
        ax.step(bins_e, np.r_[h, h[-1]], where="post",
                color=color, linestyle="--", linewidth=1.5, zorder=2,
                label=label)

    h_sig, _ = np.histogram(signal.e_p,     bins=bins_e)
    h_bkg, _ = np.histogram(background.e_p, bins=bins_e)
    diff = h_sig - h_bkg
    ax.step(bins_e, np.r_[diff, diff[-1]], where="post",
            color="black", linestyle="-", linewidth=1.5, label="Difference")

    ax.set_xlim(bins_e[0], bins_e[-1])
    ax.set_ylim(bottom=0, top=histmax * 1.25)
    _decorate_ax(ax, r"$E_p$ (MeV)", r"Entries")
    fig.tight_layout()

    # ── 2. Energy: MC Li9 vs He8 ──────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 6))
    histmax = 0.0

    for mask, color, label in [
        (mc.element == "Li9", "#ff6464", r"$^9$Li"),
        (mc.element == "He8", "#84d040", r"$^8$He"),
    ]:
        h, _ = np.histogram(mc.e_p[mask], bins=bins_e)
        histmax = max(histmax, h.max())
        ax.fill_between(bins_e, np.r_[h, h[-1]], step="post",
                        color=color, alpha=0.075, zorder=1)
        ax.step(bins_e, np.r_[h, h[-1]], where="post",
                color=color, linestyle="-", linewidth=1.5, zorder=2,
                label=label)

    ax.set_xlim(bins_e[0], bins_e[-1])
    ax.set_ylim(bottom=0, top=histmax * 1.25)
    _decorate_ax(ax, r"$E_p$ (MeV)", r"Entries")
    fig.tight_layout()

    # ── 3. dt: signal vs background ───────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 6))
    histmax = 0.0

    for data, color, label in [
        (signal.dt_mu2p,     "#648fff", "Cosmogenic enriched region"),
        (background.dt_mu2p, "#ff6464", "Cosmogenic depleted region"),
    ]:
        h, _ = np.histogram(data, bins=bins_dt)
        histmax = max(histmax, h.max())
        ax.fill_between(bins_dt, np.r_[h, h[-1]], step="post",
                        color=color, alpha=0.075, zorder=1)
        ax.step(bins_dt, np.r_[h, h[-1]], where="post",
                color=color, linestyle="--", linewidth=1.5, zorder=2,
                label=label)

    ax.set_xlim(bins_dt[0], bins_dt[-1])
    ax.set_ylim(bottom=0, top=histmax * 1.25)
    _decorate_ax(ax, r"$\Delta t$ (s)", r"Entries")
    fig.tight_layout()

    plt.show()


# ══════════════════════════════════════════════════════════════════════════════
# I/O helpers
# ══════════════════════════════════════════════════════════════════════════════

def load_data(path: str) -> tuple[CosmogenicData, CosmogenicData]:
    """Load signal and background trees from a data ROOT file."""
    fin        = uproot.open(path)
    background = CosmogenicData.from_tree(fin["background_events"])
    signal     = CosmogenicData.from_tree(fin["signal_events"])
    return signal, background


def load_simulation(path: str) -> SimulationData:
    """Load MC cosmogenic tree from a simulation ROOT file."""
    fin = uproot.open(path)
    return SimulationData.from_tree(fin["cosmogenics"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Cosmogenic isotope chi2 fit (Li9 / He8)"
    )
    parser.add_argument("--data",       type=str, required=True,
                        help="Path to data ROOT file")
    parser.add_argument("--simulation", type=str, required=True,
                        help="Path to simulation ROOT file")
    return parser.parse_args()


# ══════════════════════════════════════════════════════════════════════════════
# Entry point
# ══════════════════════════════════════════════════════════════════════════════

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

if __name__ == "__main__":
    set_latex_style()
    args = parse_args()

    # ── Load data ─────────────────────────────────────────────────────────────
    signal, background = load_data(args.data)
    mc                 = load_simulation(args.simulation)

    # ── Binning ───────────────────────────────────────────────────────────────
    bins_e_raw  = nmo_analysis_bins() # np.linspace(0.0, 12.0, 51)
    bins_dt_raw = np.linspace(0.0,  2.0, 101)
    bins_e_fit  = nmo_analysis_bins() # np.linspace(0.0, 12.0, 51)
    bins_dt_fit = np.linspace(0.0,  2.0, 101)

    # ── Pre-fit overview plots ─────────────────────────────────────────────────
    plot_raw_spectra(signal, background, mc, bins_e_raw, bins_dt_raw)

    # ── Build templates ───────────────────────────────────────────────────────
    templates = FitTemplates.build(signal, background, mc,
                                   bins_e_fit, bins_dt_fit)

    # ── Run fit ───────────────────────────────────────────────────────────────
    fitter = CosmogenicFitter(
        templates,
        fixed={},                                        # all free
        # fixed={"tau_Li9": 0.257, "tau_He8": 0.172},   # PDG-fixed lifetimes
    )
    result = fitter.fit()
    print(result.summary())

    # ── Diagnostic plots ──────────────────────────────────────────────────────
    plotter = FitPlotter(result, templates, bins_e_fit, bins_dt_fit)
    plotter.plot_all()