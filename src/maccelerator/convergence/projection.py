"""Convergence criteria that project onto a grid."""

import logging

import numpy as np
from scipy import interpolate

from .base import ConvergenceChecker
from .base import distribution_norm_tvd


__author__ = 'harrigan'

# Boltzmann constant in md units
KB = 0.0083145

log = logging.getLogger(__name__)


class PopulationProjectionTVD(ConvergenceChecker):
    """Project equilibrium populations onto a grid and compute TVD.

    :param modeller: Contains the estimated model
    :param grid: Grid onto which we project
    :param potentialfunc: Analytical (2d) potential function from which we
                          compute reference probabilities
    :param temp: Temperature at which we compute reference probabilities
    """

    def __init__(self, tolerance, grid, potentialfunc, temp):
        super().__init__(tolerance)

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
        log.debug("TVD Error: %g\t Threshold: %g", errorval, self.tolerance)

        # Return whether converged
        return errorval < self.tolerance

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
                   aspect='auto', origin='lower')
        top.scatter(sstate[:, 0], sstate[:, 1], s=100, c='w', linewidths=2,
                    zorder=10)

        # TODO: Make this into one plot
        bot.set_title("Populations (Reference)")
        bot.imshow(self.ref.T, interpolation='nearest', extent=bounds,
                   aspect='auto', origin='lower')

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
        log.debug('Making grid with bounds %.2f %.2f %.2f %.2f', *volume.bounds)

        grid_width = max(xmax - xmin, ymax - ymin) / resolution
        log.debug("Gridwithd %f", grid_width)

        grid = np.mgrid[xmin: xmax: grid_width, ymin: ymax: grid_width]
        return grid

