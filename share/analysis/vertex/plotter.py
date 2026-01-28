import argparse

import matplotlib.pyplot as plt
import numpy as np
import uproot

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ibd-analysis", type=str, nargs="+", help="IBD analysis filepath")
    parser.add_argument("--cosmo-shape-analysis", type=str, nargs="+", help="Cosmo shape analysis filepath")
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

def ibd_analysis_plot(filepath: str):
    file = uproot.open(filepath)
    tree = file["events"]

    branches = [
        "posx_p", "posy_p", "posz_p", "sec_p", "nsec_p", "e_p",
        "posx_d", "posy_d", "posz_d", "sec_d", "nsec_d", "e_d"
    ]
    data = tree.arrays(branches, library="np")

    print(f"Loaded {len(data['posx_p'])} events in {filepath}")

    # Histogram
    bins = nmo_analysis_bins()
    hist, edges = np.histogram(data["e_p"], bins=bins)
    err = np.sqrt(hist)

    centers = 0.5 * (edges[1:] + edges[:-1])
    widths = edges[1:] - edges[:-1]

    fig, ax = plt.subplots(figsize=(7, 6))

    ax.fill_between(edges, np.r_[hist, hist[-1]], step="post", color="#eff3ff", zorder=1)
    ax.errorbar(
        centers, 
        hist, 
        yerr=err, 
        xerr=widths / 2, 
        fmt="o", 
        color="#648fff", 
        markersize=4.5, 
        # markeredgewidth=1.2,
        linewidth=1.2,
        elinewidth=1.2, 
        capsize=2, 
        zorder=3
    )

    ax.set_xlabel("Prompt Energy (MeV)")
    ax.set_ylabel("Entries")

    ax.tick_params(direction="in", which="both", top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    ax.set_xlim(0.0, 12.5)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    args = parse_args()
    if args.ibd_analysis:
        for filepath in args.ibd_analysis:
            ibd_analysis_plot(filepath)
    if args.cosmo_shape_analysis:
        for filepath in args.cosmo_shape_analysis:
            pass

