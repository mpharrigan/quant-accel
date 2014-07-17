"""
Created on Mar 24, 2014

@author: harrigan

From an already built msm, check convergence and write to a file
"""

import numpy as np
import logging as log
from scipy import interpolate
from matplotlib import pyplot as plt

# Boltzmann constant in md units
KB = 0.0083145


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

    def check_convergence(self, params):
        """Check convergence

        :param params: Parameters
        :returns: bool -- whether it's converged.
        """
        raise NotImplementedError

    def plot(self, axs, sstate):
        """Plot on given axes.

        :param axs: Axes
        :param sstate: Starting states (for visualization)
        """
        raise NotImplementedError

    def plot_and_save(self, params, sstate):
        """Plot and save a visualization of the convergence check

        This is used if we *aren't* using a HybridConvergenceChecker, and
        it will just call the plot method with two axes

        :param params: Simulation parameters
        :param sstate: Starting states (to assist in visualization)
        """
        fig, axs = plt.subplots(nrows=1, ncols=self.n_plots, squeeze=True)
        self.plot(axs, sstate)

        fig.set_size_inches(7 * self.n_plots, 5)
        fig.savefig("{}.png".format(params.plot_fn))

    @property
    def n_plots(self):
        """Number of plots a particular checker generates."""
        raise NotImplementedError


class HybridConvergenceChecker(ConvergenceChecker):
    """Combine two convergence checkers with logical and"""

    def __init__(self, *checkers):
        super().__init__()
        self.checkers = checkers
        self.n_checkers = len(checkers)

    def check_convergence(self, params):
        """Check the convergence of each checker combining with logical and

        :param params: Simulation parameters
        :returns converged: bool -- whether all are converged.
        """
        converged = True
        for checker in self.checkers:
            converged = (converged and checker.check_convergence(params))

        return converged

    def plot_and_save(self, params, sstate):
        """Plot a visualization for each subchecker, stacked horizontally

        :param params: Simulation parameters
        :param sstate: Starting states
        """
        fig, axs = plt.subplots(nrows=self.n_plots, ncols=self.n_checkers,
                                squeeze=False)
        for i, checker in enumerate(self.checkers):
            checker.plot(axs[:, i], sstate)

        fig.set_size_inches(7 * self.n_checkers, 5 * self.n_plots)
        fig.suptitle(params.pretty_desc)
        fig.savefig("{}.png".format(params.plot_fn))

    @property
    def n_plots(self):
        """Find the maximum number of plots for any individual checker.

        Use this as the number of rows; Each checker will get a column of
        axes"""
        return np.max([checker.n_plots for checker in self.checkers])


class PopulationProjectionTVD(ConvergenceChecker):
    """Project equilibrium populations onto a grid and compute TVD.

    :param modeller: Contains the estimated model
    :param grid: Grid onto which we project
    :param potentialfunc: Analytical (2d) potential function from which we
                          compute reference probabilities
    :param temp: Temperature at which we compute reference probabilities
    """

    def __init__(self, modeller, grid, potentialfunc, temp):
        super().__init__()

        self.modeller = modeller
        self.grid = grid
        self.potentialfunc = potentialfunc
        self.temp = temp

        self.distribution_norm = distribution_norm_tvd

        self.ref = None
        self.est = None

    def check_convergence(self, params):
        """Use the latest model to check convergence

        :param params: Contains the threshold
        """

        xx, yy = self.grid
        n_states = self.modeller.msm.transmat_.shape[0]
        mapping = self.modeller.msm.mapping_
        back_map = []
        # Construct back-mapping to trim centroids
        for k, v in mapping.items():
            if v >= 0:
                back_map += [k]
        back_map = np.asarray(back_map)
        centroids = self.modeller.clusterer.cluster_centers_[back_map]

        # A little debug info and asserts
        log.debug(
            "Checking convergence. %d states (%d centroids, %d untrimmed)",
            n_states, len(centroids),
            len(self.modeller.clusterer.cluster_centers_))
        assert n_states == len(centroids)

        # Make sure our arrays have enough dimensions
        if len(centroids.shape) == 1:
            centroids = centroids[np.newaxis, :]
        known_points = np.vstack((centroids[:, 0], centroids[:, 1])).T

        # Figure out interpolation method
        if n_states >= 4:
            method = 'cubic'
        else:
            method = 'nearest'

        # Interpolate
        est = interpolate.griddata(known_points, self.modeller.msm.populations_,
                                   (xx, yy), method=method, fill_value=0.0)
        est = est.clip(min=0.0)
        est /= np.sum(est)
        self.est = est

        # Calculate actual result
        calc_eq = self.potentialfunc(xx, yy)
        calc_eq = np.exp(-calc_eq / (self.temp * KB))
        calc_eq /= np.sum(calc_eq)
        self.ref = calc_eq

        # Calculate error
        errorval = self.distribution_norm(calc_eq, est)
        log.debug("TVD Error: %g\t Threshold: %g", errorval, params.threshold)

        # Return whether converged
        return errorval < params.threshold

    def plot(self, axs, sstate):
        """Plot a projection onto a grid.

        :param axs: At least two axes on which we will plot
        :param sstate: Starting states to overlay
        """
        # Get axes
        top, bot = axs[0:2]

        # Set up Grid
        xx, yy = self.grid
        bounds = (xx.min(), xx.max(), yy.min(), yy.max())

        top.set_title('Populations (Estimated)')
        top.imshow(self.est.T, interpolation='nearest', extent=bounds,
                   aspect='auto',
                   origin='lower')
        top.scatter(sstate[:, 0], sstate[:, 1], s=100, c='w', linewidths=2,
                    zorder=10)

        # TODO: Make this into one plot
        bot.set_title("Populations (Reference)")
        bot.imshow(self.ref.T, interpolation='nearest', extent=bounds,
                   aspect='auto', origin='lower')

    @property
    def n_plots(self):
        return 2


class PopulationCentroidTVD(ConvergenceChecker):
    """Compare equilibrium populations of states to a reference MSM."""

    def __init__(self, modeller, centers, ref_msm):
        super().__init__()
        self.modeller = modeller
        self.centers = centers
        self.ref_msm = ref_msm

        self.distribution_norm = distribution_norm_tvd
        self.cmap = plt.get_cmap('RdBu')

        self.errors_over_time = []

    def check_convergence(self, params):
        """Check convergence

        :param params: Parameters
        :returns: Boolean, whether it has converged
        """

        est = self.modeller.full_populations
        ref = self.ref_msm.populations_

        errorval = self.distribution_norm(ref, est)
        log.debug("Population Error: %g\t Threshold: %g", errorval,
                  params.threshold)

        self.errors_over_time += [errorval]

        return errorval < params.threshold

    def plot(self, axs, sstate):
        top, bot = axs[0:2]

        est = self.modeller.msm.populations_
        ref = self.ref_msm.populations_
        scale = 1e4

        top.scatter(self.centers[:, 0], self.centers[:, 1], s=scale * ref,
                    c=ref, cmap=self.cmap, linewidths=0)
        top.scatter(self.centers[:, 0], self.centers[:, 1], s=scale * est,
                    facecolors='none', edgecolors='k', linewidths=2)
        top.set_title('Populations')

        # Bottom (over time)
        bot.plot(self.errors_over_time, 'o-')
        bot.set_ylim(0, 1)
        bot.set_xlabel('Time')

    @property
    def n_plots(self):
        return 2


class EigenvecCentroid(ConvergenceChecker):
    """Compare the first eigenvector to that of a reference MSM."""

    def __init__(self, modeller, centers, ref_msm):
        super().__init__()
        self.modeller = modeller
        self.centers = centers
        self.ref_msm = ref_msm

    def check_convergence(self, params):
        return True

    def plot(self, axs, sstate):
        top, bot = axs[0:2]

        top.scatter(self.centers[:, 0], self.centers[:, 1], s=200)
        top.set_title('Eigenvec')

        x = np.linspace(0, 1)
        bot.plot(x, 1 - x)
        bot.set_xlabel('Time')

    @property
    def n_plots(self):
        return 2


class Volume(object):
    """An object representing a rectangular area

    Use this to keep track of how big of a grid we need to make to
    project onto without discarding any info.
    """

    def __init__(self, centroidx=None, centroidy=None):
        if centroidx is not None and centroidy is not None:
            self.xmin = np.min(centroidx)
            self.xmax = np.max(centroidx)
            self.ymin = np.min(centroidy)
            self.ymax = np.max(centroidy)
        else:
            self.xmin = 0.0
            self.xmax = 0.0
            self.ymin = 0.3
            self.ymax = 0.3

    @property
    def volume(self):
        """The volume."""
        return (self.xmax - self.xmin) * (self.ymax - self.ymin)

    @property
    def bounds(self):
        """Bounds of this volume."""
        return (self.xmin, self.xmax, self.ymin, self.ymax)

    def union(self, other_v):
        """Compare to other volume and turn me into something that contains
        both.
        """
        self.xmin = min(self.xmin, other_v.xmin)
        self.ymin = min(self.ymin, other_v.ymin)
        self.xmax = max(self.xmax, other_v.xmax)
        self.ymax = max(self.ymax, other_v.ymax)

    @staticmethod
    def get_grid(volume, resolution):
        """Grid up a space

        :param volume: Info on the bounds
        :param resolution: how fine the grid should be
        """
        (xmin, xmax, ymin, ymax) = volume.bounds
        log.debug('Making grid with bounds %.2f %.2f %.2f %.2f',
                  *volume.bounds)

        grid_width = max(xmax - xmin, ymax - ymin) / resolution
        log.debug("Gridwithd %f", grid_width)

        grid = np.mgrid[xmin: xmax: grid_width, ymin: ymax: grid_width]
        return grid


    #TODO: Put these where they belong
    MULLERTHRESH = 0.6
    NTL9THRESH = 0.4

