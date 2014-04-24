"""Helper functions for plotting and visualization from ipython notebook.
"""

__author__ = 'harrigan'

from matplotlib import pyplot as plt
import numpy as np
from quantaccel import speedcalc as sc

import logging

log = logging.getLogger()


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


def _col_size(n_col, min_width=1, max_width=24, tau=4):
    """Provide slight expansion for larger columns."""
    # TODO: exponential isn't the best functional form here
    delta = max_width - min_width
    return max_width - (delta * np.exp(-n_col / tau))


def histograms(plotvs, whole_range, log_plot, htype='kde', shoop=False):
    """Plot histograms. Creates a grid of plots based on the structure of
    plotvs.

    :param shoop: roll up to one row with multiple hists per plot

    """

    if shoop:
        n_row = 1
        row_size = 5
    else:
        n_row = max([len(set(v[:, 0])) for v in plotvs.values()])
        row_size = 3 * n_row

    n_col = len(plotvs.keys())
    col_size = _col_size(n_col)
    fig, axes = plt.subplots(n_row, n_col, figsize=(col_size, row_size))
    for i, (k, v) in enumerate(plotvs.items()):
        for j, (hist, bin_centers, xval) in enumerate(
                sc.histogram(v, htype=htype, whole_range=whole_range,
                             log=log_plot)):

            # Subplots doesn't return a consistently sized array
            if n_row == n_col == 1:
                ax = axes
            elif n_row == 1:
                ax = axes[i]
            elif n_col == 1:
                ax = axes[j]
            else:
                ax = axes[j, i]

            # Try to int-ify xval
            try:
                xval = int(xval)
            except ValueError:
                pass

            # Make histogram
            fmt = '-' if htype == 'kde' else 'o-'
            ax.plot(bin_centers, hist, fmt, label=xval)
            hist_xlabel = plotvs.ylabel
            if log_plot:
                hist_xlabel = "log({})".format(hist_xlabel)
            ax.set_xlabel(hist_xlabel)

            if shoop:
                # Add legend and abbreviated title
                handles, labels = ax.get_legend_handles_labels()
                if len(labels) > 5:
                    pass
                else:
                    ax.legend(loc='best', title=plotvs.xlabel)
                ax.set_title("%s: %s" % (plotvs.labellabel, k))
            else:
                # Add complete title
                ax.set_title(
                    "%s: %s %s: %s" % (
                        plotvs.labellabel, k, plotvs.xlabel, xval))


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


def muller_non_adaptive(muller):
    rounds_vs_spt = muller.plot_vs('pop', 'rounds', xaxis='spt', label='tpr')

    # Find SPT's of interest
    maxspts = dict()
    for k, v in rounds_vs_spt.items():
        maxspts[k] = np.max(v[:, 0])
        log.debug("tpr", k, "has max spt =", np.max(v[:, 0]))

    def where_func(r):
        # Why not hack in the following while we're at it?
        r.params['adaptive'] = False
        return r.params['spt'] == maxspts[r.params['tpr']]

    # Get our dicts
    time_vs_spt_na = muller.plot_vs('sub', 'time', xaxis='spt', label='tpr',
                                    where=where_func)
    time_vs_tpr_na = muller.plot_vs('sub', 'time', xaxis='tpr',
                                    label='adaptive', where=where_func)

    log.info('; '.join(["Tpr {}: {}".format(k, v.shape[0]) for (k, v) in
                        time_vs_spt_na.items()]))
    return time_vs_spt_na, time_vs_tpr_na


def na_plot(time_vs_tpr_na):
    # TODO: Incorperate this with other plot function?
    #            - At least do most of the plotting with an auxillary function
    #              called by both.
    #            - Then add on 'perfect scaling' lines

    fig, ax = plt.subplots()
    v = time_vs_tpr_na[False]
    naw = sc.avg_and_errorbar(v)

    # Scatter
    ax.scatter(v[:, 0], v[:, 1], alpha=0.5, s=40)

    # Avgs
    ax.errorbar(naw[:, 0], naw[:, 1], yerr=np.transpose(naw[:, 2:4]), fmt='ro-',
                elinewidth=2.0, capsize=5.0, capthick=2.0, label='Observed')

    ax.set_xscale('log');
    ax.set_yscale('log')
    ax.set_xlim((1e-1, 1e4))
    ax.set_xlabel('tpr');
    ax.set_ylabel('Time')
    ax.set_title('Non-adaptive scaling')
    ylim = ax.get_ylim()

    # Perfect scaling
    from_beg = naw[0, 1] / naw[:, 0]
    from_end = (naw[-1, 1] * naw[-1, 0]) / naw[:, 0]
    ax.plot(naw[:, 0], from_beg, 'g--', label='Perfect scaling')
    ax.plot(naw[:, 0], from_end, 'g--')

    ax.set_ylim(ylim)
    ax.legend(loc='best')

    return naw


def speed_figs(plotvs, naw):
    sfig, sax = plt.subplots()
    tfig, tax = plt.subplots()
    for fit in plotvs.enum_fit():
        cfig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))

        # Scatter points
        ax1.scatter(fit.x, fit.y, s=40, alpha=0.5)

        # Plot Rounds = 1
        ax1.axhline(1.0, c='grey')

        # Plot time on right and on top fig
        dtavg = sc.avg_and_errorbar(np.array([fit.x, fit.dat_time]).T)
        ax2.scatter(fit.x, fit.dat_time, s=40, alpha=0.5)
        ax2.errorbar(dtavg[:, 0], dtavg[:, 1], yerr=np.transpose(dtavg[:, 2:4]),
                     fmt='ro-', elinewidth=1.0, capsize=1.0, capthick=1.0)

        # References
        ref_cur = naw[np.where(naw[:, 0] == fit.key)[0], 1]
        ref_fir = naw[np.where(naw[:, 0] == 1)[0], 1]

        ax2.axhline(ref_cur, c='b')
        ax2.axhline(ref_fir, c='r')

        # Calc speedup
        dt_cur = sc.avg_and_errorbar(
            np.array([fit.x, ref_cur / fit.dat_time]).T)
        dt_ref = sc.avg_and_errorbar(
            np.array([fit.x, ref_fir / fit.dat_time]).T)

        sax.errorbar(dt_cur[:, 0], dt_cur[:, 1],
                     yerr=np.transpose(dt_cur[:, 2:4]),
                     elinewidth=1.0, capsize=1.0, capthick=1.0, label=fit.key)

        tax.errorbar(dt_ref[:, 0], dt_ref[:, 1],
                     yerr=np.transpose(dt_ref[:, 2:4]),
                     elinewidth=1.0, capsize=1.0, capthick=1.0, label=fit.key)

        # Format
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel(plotvs.xlabel)
        ax1.set_ylabel(plotvs.ylabel)
        ax1.set_title('%s: %s' % (plotvs.labellabel, fit.key))

        ax2.set_xlabel(plotvs.xlabel)
        ax2.set_ylabel('Time')
        ax2.set_xscale('log')
        ax2.set_yscale('log')

    sax.set_xlabel(plotvs.xlabel)
    sax.set_ylabel('Speedup')
    sax.set_xscale('log')
    sax.set_yscale('log')
    sax.legend(loc='best')
    sax.axhline(1.0, c='grey')

    tax.set_xlabel(plotvs.xlabel)
    tax.set_ylabel('Speedup to one long traj')
    tax.set_xscale('log')
    tax.set_yscale('log')
    tax.legend(loc='best')
    tax.axhline(1.0, c='grey')