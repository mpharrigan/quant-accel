"""Base classes for convergence checker objects."""

import numpy as np
from matplotlib import pyplot as plt

__author__ = 'harrigan'


def distribution_norm_tvd(p, q):
    """Total variation distance.
    :param p: Exact
    :param q: Estimate
    """
    q /= np.sum(q)
    p /= np.sum(p)

    res = 0.5 * np.sum(np.abs(p - q))
    return res


class ConvergenceChecker:
    """Base class for checking convergence and visualizing."""

    def __init__(self, tolerance):
        self.errors_over_time = []
        self.tolerance = tolerance

    def check_convergence(self, params):
        """Check convergence

        :param params: Parameters
        :returns: bool -- whether it's converged.
        """
        raise NotImplementedError

    def fallback(self, params):
        """Called if a model could not be built."""
        self.errors_over_time += [None]

    def plot(self, axs, sstate):
        """Plot on given axes.

        :param axs: Axes
        :param sstate: Starting states (for visualization)
        """
        raise NotImplementedError

    def fallback_plot(self, axs, sstate):
        """Plot something if we couldn't build a model."""
        pass

    def plot_and_save(self, params, sstate, plot_fn, fallback=False):
        """Plot and save a visualization of the convergence check

        This is used if we *aren't* using a HybridConvergenceChecker, and
        it will just call the plot method with two axes

        :param params: Simulation parameters
        :param sstate: Starting states (to assist in visualization)
        """
        if self.n_plots < 1:
            return

        fig, axs = plt.subplots(nrows=1, ncols=self.n_plots, squeeze=True)

        if fallback:
            self.fallback_plot(axs, sstate)
        else:
            self.plot(axs, sstate)

        fig.set_size_inches(7 * self.n_plots, 5)
        fig.savefig(plot_fn)
        plt.close(fig)

    @property
    def n_plots(self):
        """Number of plots a particular checker generates."""
        raise NotImplementedError

    @property
    def plotfn(self):
        return "plot-{round_i:04d}.png"


class HybridConvergenceChecker(ConvergenceChecker):
    """Combine two convergence checkers with logical and"""

    def __init__(self, *checkers):
        super().__init__(tolerance=-1)
        self.checkers = checkers
        self.n_checkers = len(checkers)

    def check_convergence(self, params):
        """Check the convergence of each checker combining with logical and

        :param params: Simulation parameters
        :returns converged: bool -- whether all are converged.
        """
        converged = []
        for checker in self.checkers:
            # Don't do this in a more clever way or short-circuit evaluation
            # Might bite you
            converged += [checker.check_convergence(params)]

        return np.all(converged)

    def plot_and_save(self, params, sstate, plot_fn, fallback=False):
        """Plot a visualization for each subchecker, stacked horizontally

        :param params: Simulation parameters
        :param sstate: Starting states
        """
        if self.n_plots < 1:
            return

        fig, axs = plt.subplots(nrows=self.n_plots, ncols=self.n_checkers,
                                squeeze=False)
        for i, checker in enumerate(self.checkers):
            if fallback:
                checker.fallback_plot(axs[:, i], sstate)
            else:
                checker.plot(axs[:, i], sstate)

        fig.set_size_inches(7 * self.n_checkers, 5 * self.n_plots)
        fig.suptitle(params.pretty_desc)
        fig.savefig(plot_fn)
        plt.close(fig)

    @property
    def n_plots(self):
        """Find the maximum number of plots for any individual checker.

        Use this as the number of rows; Each checker will get a column of
        axes"""
        return np.max([checker.n_plots for checker in self.checkers])

