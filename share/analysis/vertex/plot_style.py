from matplotlib.ticker import AutoLocator, AutoMinorLocator, MultipleLocator

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