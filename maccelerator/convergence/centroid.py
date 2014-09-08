"""Convergence criteria that depend on having a consistent state definition."""

import logging

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import Normalize
import scipy.sparse
import scipy.linalg

from .base import SubConvergenceChecker, SubConvergence
from .base import distribution_norm_tvd

log = logging.getLogger(__name__)

__author__ = 'harrigan'
CMAP = plt.get_cmap('RdBu')


def scatter_eigenvector(ax, centers, est, ref, scale=2e3):
    """Scatter point sizes are the absolute value of eigenvector values

    Sign is indicated by color.

    :param ax: The axes object
    :param est: Estimated eigenvector values
    :param ref: Reference eigenvector values
    :param scale: Apply this multiplicative factor to the size of the points
    """
    rpos = np.where(ref >= 0)[0]
    rneg = np.where(ref < 0)[0]
    epos = np.where(est >= 0)[0]
    eneg = np.where(est < 0)[0]

    # Reference
    ax.scatter(centers[rpos, 0], centers[rpos, 1], linewidths=0,
               s=scale * ref[rpos], c='coral')
    ax.scatter(centers[rneg, 0], centers[rneg, 1], linewidths=0,
               s=-scale * ref[rneg], c='skyblue')

    # Estimate
    ax.scatter(centers[epos, 0], centers[epos, 1], linewidths=2,
               facecolors='none', edgecolors='r', s=scale * est[epos])
    ax.scatter(centers[eneg, 0], centers[eneg, 1], linewidths=2,
               facecolors='none', edgecolors='b', s=-scale * est[eneg])


class CentroidConvergenceChecker(SubConvergenceChecker):
    """For when we have consistent state decompositions."""

    def __init__(self, tolerance, centers, ref_msm):
        super().__init__(tolerance)
        self.centers = centers
        self.ref_msm = ref_msm


class PopulationCentroidTVDConvergence(SubConvergence):
    """Convergence object for PopulationCentroidTVD."""

    def __init__(self, converged, errors_over_time, est, ref, centers):
        super().__init__(converged, errors_over_time)
        self.est = est
        self.ref = ref
        self.centers = centers

    def plot(self, axs, sstate):
        """Scatter points where size shows populations

        Reference are shown as fill-only and estimate is edge-only circles.

        :param axs: Axes objects
        :param sstate: Starting states for visualization
        """
        top, bot = axs[0:2]

        est = self.est
        ref = self.ref
        scale = 1e4

        top.scatter(self.centers[:, 0], self.centers[:, 1], s=scale * ref,
                    c=ref, cmap=CMAP, linewidths=0, norm=Normalize(vmin=0))
        top.scatter(self.centers[:, 0], self.centers[:, 1], s=scale * est,
                    facecolors='none', edgecolors='k', linewidths=2)
        top.set_title('Populations')

        self._plot_bottom(bot)

    @property
    def n_plots(self):
        return 2


class PopulationCentroidTVD(CentroidConvergenceChecker):
    """Compare equilibrium populations of states to a reference MSM."""

    def __init__(self, tolerance, centers, ref_msm):
        super().__init__(tolerance, centers, ref_msm)
        self.distribution_norm = distribution_norm_tvd

    def check_convergence(self, model, params):
        """Check convergence

        :param model: The model for which we will check the convergence
        :param params: Parameters
        """

        est = model.populations
        ref = self.ref_msm.populations_

        errorval = self.distribution_norm(ref, est)
        log.debug("Population Error:\t%g\tThreshold: %g", errorval,
                  self.tolerance)

        self.errors_over_time += [errorval]
        converged = errorval < self.tolerance

        return PopulationCentroidTVDConvergence(converged,
                                                self.errors_over_time, est, ref,
                                                self.centers)


class EigenvecCentroidConvergence(SubConvergence):
    """Return object corresponding to EigenvecCentroid."""

    def __init__(self, converged, errors_over_time, applied_eigenv, ref_eigenv,
                 centers):
        super().__init__(converged, errors_over_time)

        self.applied_eigenv = np.squeeze(applied_eigenv)
        self.ref_eigenv = ref_eigenv
        self.centers = centers

    def plot(self, axs, sstate):
        """Scatter point sizes are the absolute value of eigenvector values

        Sign is indicated by color.

        :param axs: Axes objects
        :param sstate: Starting states for visualization
        """
        top, bot = axs[0:2]

        est = self.applied_eigenv
        ref = self.ref_eigenv

        scatter_eigenvector(top, self.centers, est, ref)
        top.set_title(r'$<\phi_1 | \hat{T} \circ \phi_1>$')

        self._plot_bottom(bot)

    @property
    def n_plots(self):
        return 2


class EigenvecCentroid(CentroidConvergenceChecker):
    """Compare the first eigenvector to that of a reference MSM."""

    def check_convergence(self, model, params):
        """Apply estimated tmat to reference eigenvector

        Dot with the reference eigenvector, normalize and compare with
        the reference eigenvalue.

        <y | yT> = <y | l | y'> = l <y|y'> ~ l

        :param model: The model for which we check convergence
        :param params: Run parameters (unused)
        """
        ref = self.ref_msm.left_eigenvectors_[:, 1]
        refs = scipy.sparse.csr_matrix(ref)
        tmat = model.tmat

        ref_val = self.ref_msm.eigenvalues_[1]

        # <y, Poy> / <y, y>
        applied_eigenv = refs.dot(tmat)
        est_val = refs.dot(applied_eigenv.T)[0, 0] / np.dot(ref, ref)

        errorval = ref_val - est_val
        self.errors_over_time += [errorval]

        log.debug("Eigenvec Error:\t%g\tThreshold: %g", errorval,
                  self.tolerance)

        converged = errorval < self.tolerance
        return EigenvecCentroidConvergence(converged, self.errors_over_time,
                                           applied_eigenv, ref, self.centers)


class EigenvecL2Convergence(SubConvergence):
    """Return object corresponding to EigenvecL2"""

    def __init__(self, converged, errors_over_time, ref_eigenv, est_eigenv,
                 centers):
        super().__init__(converged, errors_over_time)
        self.centers = centers
        self.ref_eigenv = ref_eigenv
        self.est_eigenv = est_eigenv

    def plot(self, axs, sstate):
        top, bot = axs[0:2]

        ref = self.ref_eigenv
        est = self.est_eigenv

        scatter_eigenvector(top, self.centers, est, ref)
        top.set_title(r'$||\hat{\phi}_1-\phi_1||$')

        self._plot_bottom(bot)

    @property
    def n_plots(self):
        return 2


class EigenvecL2(CentroidConvergenceChecker):
    """Try just doing L2 norm between eigenvector values."""

    def check_convergence(self, model, params):
        ref = self.ref_msm.left_eigenvectors_[:, 1]
        est = model.eigenvectors

        # Only relative sign matters, pick the best.
        diff1 = ref - est
        diff2 = ref + est
        errorval = min(np.dot(diff1, diff1), np.dot(diff2, diff2))
        self.errors_over_time += [errorval]

        log.debug("Eigenvec L2 Error: %g\t Threshold: %g", errorval,
                  self.tolerance)

        converged = errorval < self.tolerance
        return EigenvecL2Convergence(converged, self.errors_over_time, ref, est,
                                     self.centers)


class TMatFroConvergence(SubConvergence):
    """Return object corresponding to TMatFro"""

    def __init__(self, converged, errors_over_time):
        super().__init__(converged, errors_over_time)

    def plot(self, axs, sstate):
        top, bot = axs[0:2]

        top.set_title('Tmat Frobenius Norm')

        self._plot_bottom(bot)

    @property
    def n_plots(self):
        return 2


class TMatFro(CentroidConvergenceChecker):
    """Frobenius norm."""

    def check_convergence(self, model, params):
        est = model.tmat
        ref = self.ref_msm.transmat_

        errorval = scipy.linalg.norm(est - ref, ord='fro')
        self.errors_over_time += [errorval]

        log.debug("Frobenius Norm:\t%g\tThreshold: %g", errorval,
                  self.tolerance)

        converged = errorval < self.tolerance
        return TMatFroConvergence(converged, self.errors_over_time)


