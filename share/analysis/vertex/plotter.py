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
    bins = np.array([])
    for k in range(len(bins)):
        start = edges[k]
        end = edges[k+1]
        nbin = nbins[k]
        if k == 0:
            bins = np.linspace(start, end, nbin)
        else:
            bins = np.append(bins, np.linspace(start, end, nbin))
    return bins
        

def ibd_analysis_plot(filepath : str):
    file = uproot.open(filepath)
    tree = file["events"]
    branches = [
        "posx_p", "posy_p", "posz_p", "sec_p", "nsec_p", "e_p",
        "posx_d", "posy_d", "posz_d", "sec_d", "nsec_d", "e_d"
    ]
    data = tree.arrays(branches, library="np")

    print(f"Loaded{len(data['posx_p'])} events in {filepath}")

    fig_e_p, ax_e_p = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    hist_e_p,edges_e_p = np.histogram(data["e_p"], bins=nmo_analysis_bins())
    histerr_e_p = np.sqrt(hist_e_p)
    ax_e_p.errorbar(edges_e_p[:-1], hist_e_p, yerr=histerr_e_p, fmt="o")
    plt.show()

if __name__ == "__main__":
    args = parse_args()
    if args.ibd_analysis:
        for filepath in args.ibd_analysis:
            ibd_analysis_plot(filepath)
    if args.cosmo_shape_analysis:
        for filepath in args.cosmo_shape_analysis:
            pass

