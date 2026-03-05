#!/usr/bin/env python3

import datetime
import math
import time as pytime
from typing import Tuple

import matplotlib as mpl
import numpy as np
from scipy.spatial import cKDTree

# ---------------------------
# Utility functions
# ---------------------------

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

def generate_grid_points_in_sphere(radius : float, spacing : float, maxpoints : int = None) -> np.ndarray:
    coords = np.arange(-radius, radius, spacing)
    x, y, z = np.meshgrid(coords, coords, coords, indexing="xy")
    points = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
    r2 = np.sum(points**2, axis=1)
    inside = r2 <= radius**2
    points = points[inside]
    if maxpoints is not None and len(points) > maxpoints:
        index = np.linspace(0, len(points) - 1, maxpoints).astype(int)
        points = points[index]
    return points

def sample_times_uniform(n : int, start : float, end : float, rng : np.random.Generator) -> np.ndarray:
    return rng.uniform(start, end, size=n)

def merge_intervals(start : np.ndarray, end : np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if len(start) == 0:
        return np.array([]), np.array([])
    order = np.argsort(start)
    s = start[order].copy()
    e = end[order].copy()
    out_s = [s[0]]
    out_e = [e[0]]
    for i in range(1, len(s)):
        if s[i] <= out_e[-1]:
            out_e[-1] = max(out_e[-1], e[i])
        else:
            out_s.append(s[i])
            out_e.append(e[i])
    return np.array(out_s), np.array(out_e)

def intervals_masked(time : np.ndarray, start : np.ndarray, end : np.ndarray) -> np.ndarray:
    if len(start) == 0:
        return np.zeros(len(time), dtype=bool)
    index = np.searchsorted(start, time, side="right")
    mask = np.zeros(len(time), dtype=bool)
    valid = (index >= 0)
    mask[valid] = (time[valid] <= end[index[valid]])
    return mask

# ---------------------------
# Veto implementations
# ---------------------------

def build_muon_veto_intervals(time : np.ndarray, window : float) -> Tuple[np.ndarray, np.ndarray]:
    if len(time) == 0:
        return np.array([]), np.array([])
    start = np.asarray(time)
    end = start + window
    return merge_intervals(start, end)

def apply_muon_veto_to_testpoints(time : np.ndarray, start : np.ndarray, end : np.ndarray):
    return intervals_masked(time, start, end)

def apply_neutron_veto(test_pos: np.ndarray, test_times: np.ndarray, neutron_pos: np.ndarray, neutron_times: np.ndarray, radius_m: float = 4.0, time_window_s: float = 1.2):
    n = len(test_times)
    veto_mask = np.zeros(n, dtype=bool)
    if len(neutron_times) == 0:
        return veto_mask
    tree = cKDTree(neutron_pos)
    CH = 200_000  # chunk size (tune depending on memory)
    for i0 in range(0, n, CH):
        i1 = min(n, i0 + CH)
        pts_chunk = test_pos[i0:i1]
        t_chunk = test_times[i0:i1]
        # list of arrays of nearby neutron indices
        neighbors = tree.query_ball_point(pts_chunk, r=radius_m)
        for j, neigh in enumerate(neighbors):
            if len(neigh) == 0:
                continue
            t = t_chunk[j]
            # check if any neutron in neigh satisfies 0 < t - t_neu <= time_window_s
            dt = t - neutron_times[neigh]
            if np.any((dt > 0.0) & (dt <= time_window_s)):
                veto_mask[i0 + j] = True
    return veto_mask

def apply_track_veto(test_pos: np.ndarray, test_times: np.ndarray,
                     muon_tracks: dict,
                     radius_m: float = 1.0, time_window_s: float = 0.5):
    n = len(test_times)
    veto_mask = np.zeros(n, dtype=bool)
    if len(muon_tracks.get('times', [])) == 0:
        return veto_mask

    # We'll process test points in chunks to keep memory usage moderate.
    CH = 200_000
    times = muon_tracks['times']
    r0s = muon_tracks['r0']
    us = muon_tracks['u']
    ntracks = len(times)

    for i0 in range(0, n, CH):
        i1 = min(n, i0 + CH)
        pts_chunk = test_pos[i0:i1]  # (M,3)
        t_chunk = test_times[i0:i1]  # (M,)
        # For each track, find which points in chunk are within time window
        for k in range(ntracks):
            dt = t_chunk - times[k]
            valid = (dt > 0.0) & (dt <= time_window_s)
            if not np.any(valid):
                continue
            idxs = np.nonzero(valid)[0]
            vec = pts_chunk[idxs] - r0s[k]  # (m,3)
            # perpendicular distance to infinite line: |vec - (vec.u)u|
            proj = np.einsum('ij,j->i', vec, us[k])  # dot with u
            perp = vec - np.outer(proj, us[k])
            d2 = np.sum(perp**2, axis=1)
            veto_here = d2 < (radius_m**2)
            veto_mask[i0 + idxs[veto_here]] = True
    return veto_mask

def apply_multiplicity_veto(test_times: np.ndarray,
                            singles_times: np.ndarray,
                            singles_energies: np.ndarray,
                            singles_muon_veto_mask: np.ndarray,
                            energy_window: Tuple[float,float] = (2.0, 12.0),
                            window_s: float = 1e-3,
                            tau_delay_s: float = 220e-6):
    eff_window = window_s + tau_delay_s
    # preselect singles by energy and not muon-vetoed
    sel = (~singles_muon_veto_mask) & (singles_energies >= energy_window[0]) & (singles_energies <= energy_window[1])
    times_sel = np.sort(singles_times[sel])
    if len(times_sel) == 0:
        return np.zeros(len(test_times), dtype=bool)
    # use searchsorted to count number in [t-eff, t+eff]
    left = np.searchsorted(times_sel, test_times - eff_window, side='left')
    right = np.searchsorted(times_sel, test_times + eff_window, side='right')
    counts = right - left
    return counts > 0

# ---------------------------
# Efficiency routine tying everything together
# ---------------------------

def grid_efficiency(test_pos: np.ndarray,
                    test_times: np.ndarray,
                    muon_times: np.ndarray,
                    muon_tracks: dict,
                    neutron_pos: np.ndarray,
                    neutron_times: np.ndarray,
                    singles_times: np.ndarray,
                    singles_energies: np.ndarray,
                    params: dict,
                    chunk_size:int = 200_000):
    n = len(test_times)
    assert test_pos.shape[0] == n

    # --- muon fixed-time veto intervals (merged)
    muon_intervals_s, muon_intervals_e = build_muon_veto_intervals(muon_times, params['veto_muon_window_s'])

    # apply muon veto (vectorized)
    mask_muon = apply_muon_veto_to_testpoints(test_times, muon_intervals_s, muon_intervals_e)

    # apply muon veto on singles (needed for multiplicity)
    singles_muon_mask = intervals_masked(singles_times, muon_intervals_s, muon_intervals_e)

    # multiplicity veto (vectorized over times)
    mask_mult = apply_multiplicity_veto(test_times,
                                        singles_times, singles_energies,
                                        singles_muon_mask,
                                        energy_window=params['mult_energy_window'],
                                        window_s=params['mult_window_s'],
                                        tau_delay_s=params['tau_delay_s'])

    # neutron veto (chunked KDTree approach)
    mask_neu = apply_neutron_veto(test_pos, test_times, neutron_pos, neutron_times,
                                  radius_m=params['neutron_radius_m'], time_window_s=params['neutron_time_s'])

    # track veto (chunked loop over tracks)
    mask_track = apply_track_veto(test_pos, test_times, muon_tracks,
                                  radius_m=params['track_radius_m'], time_window_s=params['track_time_s'])

    # Combined mask (any veto)
    mask_any = mask_muon | mask_mult | mask_neu | mask_track

    def eff_and_err(mask):
        Ntotal = len(mask)
        Npass = np.count_nonzero(~mask)
        eff = Npass / Ntotal
        # binomial uncertainty
        sigma = np.sqrt(eff*(1-eff)/Ntotal) if Ntotal>0 else 0.0
        return eff, sigma, Npass, Ntotal

    results = {
        'muon': eff_and_err(mask_muon),
        'multiplicity': eff_and_err(mask_mult),
        'neutron': eff_and_err(mask_neu),
        'track': eff_and_err(mask_track),
        'combined': eff_and_err(mask_any)
    }

    masks = {
        'mask_muon': mask_muon,
        'mask_mult': mask_mult,
        'mask_neu': mask_neu,
        'mask_track': mask_track,
        'mask_any': mask_any
    }

    return results, masks

# ---------------------------
# Toy-data generator for quick testing
# ---------------------------

def toy_data(n_muons=2000, n_neu=2000, n_singles=10000, run_start=0.0, run_end=3600.0, rng=None):
    """
    Create toy muon times, tracks, neutron candidates and singles.
    Positions generated uniformly inside sphere of R=16.5 m for e.g.
    """
    if rng is None:
        rng = np.random.default_rng(12345)

    # muon times uniformly in run
    muon_times = rng.uniform(run_start, run_end, size=n_muons)
    # muon tracks: random entry point on sphere and random direction
    R = 19.5  # approximate outer radius to place muon track origin (not critical)
    theta = np.arccos(rng.uniform(-1, 1, size=n_muons))
    phi = rng.uniform(0, 2*np.pi, size=n_muons)
    r0 = np.vstack((R*np.sin(theta)*np.cos(phi), R*np.sin(theta)*np.sin(phi), R*np.cos(theta))).T
    # random directions
    uvec = rng.normal(size=(n_muons, 3))
    uvec /= np.linalg.norm(uvec, axis=1)[:,None]
    muon_tracks = {'times': muon_times, 'r0': r0, 'u': uvec}

    # neutron candidates times and positions (a bit delayed after some muons)
    neu_times = muon_times[:min(len(muon_times), n_neu)] + rng.uniform(1e-5, 1e-2, size=min(len(muon_times), n_neu))
    # random positions inside sphere of radius 16.5 m
    Rfv = 16.5
    u = rng.normal(size=(len(neu_times), 3))
    u /= np.linalg.norm(u, axis=1)[:,None]
    r_neu = rng.uniform(0, 1, size=(len(neu_times),1))**(1/3) * Rfv * u

    # singles times and energies
    singles_times = rng.uniform(run_start, run_end, size=n_singles)
    singles_energies = rng.uniform(0, 20, size=n_singles)

    return muon_times, muon_tracks, r_neu, neu_times, singles_times, singles_energies

# ---------------------------
# Example main: build grid, run efficiency
# ---------------------------

if __name__ == "__main__":

    cd_sphere = 20050.0
    acrylic_radius = 17700.0
    fiducial_radius = 16500.0
    spacing = 100.0
    maxpoints = None
    seed = 12345

    rng = np.random.default_rng(seed)
    print("Generating grid...")
    t0 = pytime.time()
    points = generate_grid_points_in_sphere(fiducial_radius, spacing, maxpoints=maxpoints)
    npoints = len(points)
    t1 = pytime.time(); print(f"  -> {npoints} points generated in {t1-t0:.2f}s")

    start, end = 0.0, 3600.0
    time = sample_times_uniform(npoints, start, end, rng)

    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    from matplotlib.ticker import AutoMinorLocator

    set_latex_style()

    fig, ax = plt.subplots(1, 3, figsize=(16, 5))

    xy_slide_mask = np.abs(points[:, 2]) < 50.0
    xy_slide_points = points[xy_slide_mask] 
    ax[0].scatter(xy_slide_points[:, 0] / 1000.0, xy_slide_points[:, 1] / 1000.0, color="tab:blue", alpha=0.1, s=1)
    fv_circle = Circle((0, 0), fiducial_radius / 1000.0, edgecolor="red", facecolor="none", linestyle="--", linewidth=2, label="Fiducial Volume")
    ax[0].add_patch(fv_circle)
    acrylic_circle = Circle((0, 0), acrylic_radius / 1000.0, edgecolor="black", facecolor="none", linestyle="--", linewidth=2, label="Acrylic Sphere")
    ax[0].add_patch(acrylic_circle)
    
    ax[0].set_title(rf"XY ($|z| < 50$ mm), {len(xy_slide_points)}")
    ax[0].set_xlabel(r"$x$ (m)")
    ax[0].set_ylabel(r"$y$ (m)")
    # ax[0].legend()
    ax[0].minorticks_on()
    ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))
    ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))
    ax[0].tick_params(direction="in", which="both", top=True, right=True)
    ax[0].set_xlim(left=-cd_sphere / 1000.0, right=cd_sphere / 1000.0)
    ax[0].set_ylim(bottom=-cd_sphere / 1000.0, top=cd_sphere / 1000.0)
    ax[0].set_xscale("linear")
    ax[0].set_yscale("linear")

    xz_slide_mask = np.abs(points[:, 1]) < 50.0
    xz_slide_points = points[xz_slide_mask] 
    ax[1].scatter(xz_slide_points[:, 0] / 1000.0, xz_slide_points[:, 2] / 1000.0, color="tab:green", alpha=0.1, s=1)
    fv_circle = Circle((0, 0), fiducial_radius / 1000.0, edgecolor="red", facecolor="none", linestyle="--", linewidth=2, label="Fiducial Volume")
    ax[1].add_patch(fv_circle)
    acrylic_circle = Circle((0, 0), acrylic_radius / 1000.0, edgecolor="black", facecolor="none", linestyle="--", linewidth=2, label="Acrylic Sphere")
    ax[1].add_patch(acrylic_circle)

    ax[1].set_title(rf"XZ ($|y| < 50$ mm), {len(xz_slide_points)}")
    ax[1].set_xlabel(r"$x$ (m)")
    ax[1].set_ylabel(r"$z$ (m)")
    # ax[1].legend()
    ax[1].minorticks_on()
    ax[1].xaxis.set_minor_locator(AutoMinorLocator(5))
    ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))
    ax[1].tick_params(direction="in", which="both", top=True, right=True)
    ax[1].set_xlim(left=-cd_sphere / 1000.0, right=cd_sphere / 1000.0)
    ax[1].set_ylim(bottom=-cd_sphere / 1000.0, top=cd_sphere / 1000.0)
    ax[1].set_xscale("linear")
    ax[1].set_yscale("linear")

    yz_slide_mask = np.abs(points[:, 0]) < 50.0
    yz_slide_points = points[yz_slide_mask] 
    ax[2].scatter(yz_slide_points[:, 1] / 1000.0, yz_slide_points[:, 2] / 1000.0, color="tab:red", alpha=0.1, s=1)
    fv_circle = Circle((0, 0), fiducial_radius / 1000.0, edgecolor="red", facecolor="none", linestyle="--", linewidth=2, label="Fiducial Volume")
    ax[2].add_patch(fv_circle)
    acrylic_circle = Circle((0, 0), acrylic_radius / 1000.0, edgecolor="black", facecolor="none", linestyle="--", linewidth=2, label="Acrylic Sphere")
    ax[2].add_patch(acrylic_circle)

    ax[2].set_title(rf"YZ ($|x| < 50$ mm), {len(yz_slide_points)}")
    ax[2].set_xlabel(r"$y$ (m)")
    ax[2].set_ylabel(r"$z$ (m)")
    # ax[2].legend()
    ax[2].minorticks_on()
    ax[2].xaxis.set_minor_locator(AutoMinorLocator(5))
    ax[2].yaxis.set_minor_locator(AutoMinorLocator(5))
    ax[2].tick_params(direction="in", which="both", top=True, right=True)
    ax[2].set_xlim(left=-cd_sphere / 1000.0, right=cd_sphere / 1000.0)
    ax[2].set_ylim(bottom=-cd_sphere / 1000.0, top=cd_sphere / 1000.0)
    ax[2].set_xscale("linear")
    ax[2].set_yscale("linear")

    fig.tight_layout()
    fig.show()

    plt.show()

    # create toy veto event data
    # muon_times, muon_tracks, neu_pos, neu_times, singles_times, singles_energies = toy_data(
    #     n_muons=2000, n_neu=2000, n_singles=10000, run_start=run_start, run_end=run_end, rng=rng)

    # parameters for vetoes (tune to your configuration)
    # params = {
    #     'veto_muon_window_s': 5e-3,         # fixed muon veto (5 ms)
    #     'mult_window_s': 1e-3,              # multiplicity ±1 ms
    #     'mult_energy_window': (2.0, 12.0),  # energy window for multiplicity
    #     'tau_delay_s': 220e-6,              # effective delay to account delayed candidate
    #     'neutron_radius_m': 4.0,            # neutron veto radius
    #     'neutron_time_s': 1.2,              # neutron veto duration
    #     'track_radius_m': 1.0,              # track veto radius
    #     'track_time_s': 0.5                 # track veto duration
    # }

    # print("Running grid efficiency calculation (this may take some time)...")
    # t0 = pytime.time()
    # results, masks = grid_efficiency(pts, times,
    #                                  muon_times, muon_tracks,
    #                                  neu_pos, neu_times,
    #                                  singles_times, singles_energies,
    #                                  params)
    
    # t1 = pytime.time()
    # print(f"Done in {t1-t0:.2f}s\n")

    # print nicely
    # print("Veto efficiencies (pass fraction and binomial uncertainty):")
    # for k in ['muon','multiplicity','neutron','track','combined']:
    #     eff, sigma, Npass, Ntot = results[k]
    #     print(f"  {k:12s}: eff = {eff*100:6.3f}%  +/- {sigma*100:.3f}%  (Npass={Npass}/{Ntot})")