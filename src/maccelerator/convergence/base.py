"""Base classes for convergence checker objects."""

import pickle

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

    def check_convergence(self, model, params):
        """Check convergence

        :param params: Parameters
        :returns: bool -- whether it's converged.
        """
        raise NotImplementedError

    @property
    def convfn(self):
        return "convergence-{round_i}"

    @property
    def plotfn(self):
        return "plot-{round_i:04d}"


class IConvergence:
    """Interface only."""

    def plot_and_save(self, params, sstate, plot_fn):
        raise NotImplementedError

    def plot(self, axs, sstate):
        raise NotImplementedError

    def save(self, conv_fn):
        fn = "{}.pickl".format(conv_fn)
        with open(fn, 'wb') as f:
            pickle.dump(self, f)

    @property
    def converged(self):
        raise NotImplementedError


class Convergence(IConvergence):
    def __init__(self, converged, errors_over_time):
        self._converged = converged
        self.errors_over_time = errors_over_time

    @property
    def converged(self):
        return self._converged

    def plot(self, axs, sstate):
        """Plot on given axes.

        :param axs: Axes
        :param sstate: Starting states (for visualization)
        """
        raise NotImplementedError

    def _plot_bottom(self, bot):
        bot.plot(self.errors_over_time, 'o-')
        bot.axhline(0, c='k')
        bot.set_xlabel('Time')
        cstring = 'Converged' if self.converged else 'Not Converged'
        ccolor = 'green' if self.converged else 'red'
        bot.set_title(cstring, color=ccolor)


    def plot_and_save(self, params, sstate, plot_fn):
        """Plot and save a visualization of the convergence check

        This is used if we *aren't* using a HybridConvergenceChecker, and
        it will just call the plot method with two axes

        :param params: Simulation parameters
        :param sstate: Starting states (to assist in visualization)
        """
        if self.n_plots < 1:
            return

        fig, axs = plt.subplots(nrows=1, ncols=self.n_plots, squeeze=True)

        self.plot(axs, sstate)

        fig.set_size_inches(7 * self.n_plots, 5)
        fig.savefig(plot_fn)
        plt.close(fig)

    @property
    def n_plots(self):
        """Number of plots a particular checker generates."""
        raise NotImplementedError


class HybridConvergenceChecker(ConvergenceChecker):
    """Combine two convergence checkers with logical and"""

    def __init__(self, *checkers):
        super().__init__(tolerance=-1)
        self.checkers = checkers
        self.n_checkers = len(checkers)

    def check_convergence(self, model, params):
        """Check the convergence of each checker combining with logical and

        :param params: Simulation parameters
        :returns converged: bool -- whether all are converged.
        """
        convergences = []
        for checker in self.checkers:
            # Don't do this in a more clever way or short-circuit evaluation
            # Might bite you
            convergences += [checker.check_convergence(model, params)]

        return HybridConvergence(convergences)


class HybridConvergence(IConvergence):
    def __init__(self, convergences):
        super().__init__()
        self.convergences = convergences
        self.n_convergences = len(convergences)

    @property
    def converged(self):
        return np.all([c.converged for c in self.convergences])

    def plot_and_save(self, params, sstate, plot_fn):
        """Plot a visualization for each subchecker, stacked horizontally

        :param params: Simulation parameters
        :param sstate: Starting states
        """
        if self.n_plots < 1:
            return

        fig, axs = plt.subplots(nrows=self.n_plots, ncols=self.n_convergences,
                                squeeze=False)
        for i, convergence in enumerate(self.convergences):
            convergence.plot(axs[:, i], sstate)

        fig.set_size_inches(7 * self.n_convergences, 5 * self.n_plots)
        fig.suptitle(params.pretty_desc)
        fig.savefig("{}.png".format(plot_fn))
        plt.close(fig)

    @property
    def n_plots(self):
        """Find the maximum number of plots for any individual checker.

        Use this as the number of rows; Each checker will get a column of
        axes"""
        return np.max([conv.n_plots for conv in self.convergences])

