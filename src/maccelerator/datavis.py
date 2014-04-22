"""Helper functions for plotting and visualization from ipython notebook.
"""

__author__ = 'harrigan'

from matplotlib import pyplot as plt
import numpy as np
from quantaccel import speedcalc as sc


def convergence_plots(resobj, runcopy, name='pop', y_label='Pop TVD'):
    """Plot convergence vs. time. Default is population TVD

    This will produce a grid of plots, one for each lt-spt combo.
    """
    fstring, tuple_gen, uniques = resobj.get_unique(name, 'lt', 'spt',
                                                    'adaptive')

    # Make it into a dict for plotting
    uniques = [(ss, ax) for ax, (sv, ss) in enumerate(uniques)]
    uniqued = dict(uniques)

    # Figure out number of plots
    n_configs = len(uniques)
    n_cols = 3
    n_rows = (n_configs + 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows))
    for r in resobj.get_runcopy(runcopy, name):
        # Plot on appropriate axes
        figname = fstring % tuple_gen(r)
        ax = axes[np.unravel_index(uniqued[figname], axes.shape)]

        ax.plot(r.errors[:, 0], r.errors[:, 1], '.-',
                label=r.params[resobj.param_compat['tpr']])

        ax.set_title(figname)
        ax.set_xlabel('wall time')
        ax.set_ylabel(y_label)
        ax.axhline(0.6, c='r')
        ax.set_xscale('log')
        ax.legend(loc='best')


def histograms(plotvs, whole_range, log, htype='kde'):
    """Plot histograms. Creates a grid of plots based on the structure of
    plotvs.

    """

    n_col = len(plotvs.keys())
    n_row = max([len(set(v[:, 0])) for v in plotvs.values()])

    fig, axes = plt.subplots(n_row, n_col,
                             figsize=((16 / 5.0) * n_col, 3 * n_row))
    for i, (k, v) in enumerate(plotvs.items()):
        for j, (hist, bin_centers, xval) in enumerate(
                sc.histogram(v, htype=htype, whole_range=whole_range, log=log)):
            ax = axes[j, i]
            fmt = '-' if htype == 'kde' else 'o-'
            ax.plot(bin_centers, hist, fmt)
            ax.set_xlabel('# of rounds')
            ax.set_title("tpr: %s spt: %s" % (k, xval))


def plot_avgs(plotvs):
    """Plot all on one axis and seperate axes."""
    fig, ax = plt.subplots(figsize=(10, 7))

    n_cols = 3
    n_configs = len(plotvs.items())
    n_rows = (n_configs + 1) // n_cols

    fig2, axes = plt.subplots(n_rows, n_cols, figsize=(12, n_rows * 4))
    axes = np.ravel(axes)

    for i, (k, v) in enumerate(plotvs.items()):
        w = sc.avg_and_errorbar(v)

        # Top plot
        ax.errorbar(w[:, 0], w[:, 1], yerr=np.transpose(w[:, 2:4]), fmt='o-',
                    elinewidth=2.0, capsize=4.0, capthick=2.0, label=k)

        # Individual plot
        axes[i].errorbar(w[:, 0], w[:, 1], yerr=np.transpose(w[:, 2:4]),
                         fmt='ro-', elinewidth=2.0, capsize=4.0, capthick=2.0,
                         label=k)

        axes[i].set_xscale('log')
        axes[i].set_yscale('log')
        axes[i].set_xlabel(plotvs.xlabel)
        axes[i].set_ylabel(plotvs.ylabel)
        axes[i].axhline(1.0, c='grey')
        axes[i].set_title('tpr = %d' % k)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(plotvs.xlabel)
    ax.set_ylabel(plotvs.ylabel)
    ax.legend(loc='best', title=plotvs.labellabel)
    ax.axhline(1.0, c='grey')
