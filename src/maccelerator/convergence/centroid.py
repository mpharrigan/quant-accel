"""Convergence criteria that depend on having a consistent state definition."""

import logging as log

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import Normalize
import scipy.sparse
import scipy.linalg

from .base import ConvergenceChecker
from .base import distribution_norm_tvd


__author__ = 'harrigan'


def scatter_eigenvector(ax, centers, est, ref, scale=2e3):
    """Scatter point sizes are the absolute value of eigenvector values

    Sign is indicated by color.
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


class CentroidConvergenceChecker(ConvergenceChecker):
    """For when we have consistent state decompositions."""

    def __init__(self, tolerance, modeller, centers, ref_msm):
        super().__init__(tolerance)

        self.modeller = modeller
        self.centers = centers
        self.ref_msm = ref_msm


class PopulationCentroidTVD(CentroidConvergenceChecker):
    """Compare equilibrium populations of states to a reference MSM."""

    def __init__(self, tolerance, modeller, centers, ref_msm):
        super().__init__(tolerance, modeller, centers, ref_msm)
        self.distribution_norm = distribution_norm_tvd
        self.cmap = plt.get_cmap('RdBu')


    def check_convergence(self, params):
        """Check convergence

        :param params: Parameters
        :returns: Boolean, whether it has converged
        """

        est = self.modeller.full_populations
        ref = self.ref_msm.populations_

        errorval = self.distribution_norm(ref, est)
        log.debug("Population Error: %g\t Threshold: %g", errorval,
                  self.tolerance)

        self.errors_over_time += [errorval]

        return errorval < self.tolerance

    def plot(self, axs, sstate):
        """Scatter points where size shows populations

        Reference are shown as fill-only and estimate is edge-only circles.
        """
        top, bot = axs[0:2]

        est = self.modeller.msm.populations_
        ref = self.ref_msm.populations_
        scale = 1e4

        top.scatter(self.centers[:, 0], self.centers[:, 1], s=scale * ref,
                    c=ref, cmap=self.cmap, linewidths=0, norm=Normalize(vmin=0))
        top.scatter(self.centers[:, 0], self.centers[:, 1], s=scale * est,
                    facecolors='none', edgecolors='k', linewidths=2)
        top.set_title('Populations')

        # Bottom (over time)
        bot.plot(self.errors_over_time, 'o-')
        bot.axhline(0, c='k')
        bot.set_xlabel('Time')

    @property
    def n_plots(self):
        return 2


class EigenvecCentroid(CentroidConvergenceChecker):
    """Compare the first eigenvector to that of a reference MSM."""

    def check_convergence(self, params):
        """Apply estimated tmat to reference eigenvector

        Dot with the reference eigenvector, normalize and compare with
        the reference eigenvalue.

        <y | yT> = <y | l | y'> = l <y|y'> ~ l
        """
        ref = self.ref_msm.eigenvectors_[:, 1]
        refs = scipy.sparse.csr_matrix(ref)
        tmat = self.modeller.full_tmat

        ref_val = self.ref_msm.eigenvalues_[1]

        # <y, Poy> / <y, y>
        self.applied_eigenv = refs.dot(tmat)
        est_val = refs.dot(self.applied_eigenv.T)[0, 0] / np.dot(ref, ref)

        errorval = ref_val - est_val
        self.errors_over_time += [errorval]

        log.debug("Eigenvec Error: %g\t Threshold: %g", errorval,
                  self.tolerance)

        return errorval < self.tolerance


    def plot(self, axs, sstate):
        """Scatter point sizes are the absolute value of eigenvector values

        Sign is indicated by color.
        """
        top, bot = axs[0:2]

        est = self.applied_eigenv.toarray().flatten()
        ref = self.ref_msm.eigenvectors_[:, 1]

        scatter_eigenvector(top, self.centers, est, ref)
        top.set_title('$<\phi_1 | \hat{T} \circ \phi_1>$')

        # Bottom (over time)
        bot.plot(self.errors_over_time, 'o-')
        bot.axhline(0, c='k')
        bot.set_xlabel('Time')

    @property
    def n_plots(self):
        return 2


class EigenvecL2(CentroidConvergenceChecker):
    """Try just doing L2 norm between eigenvector values."""

    def check_convergence(self, params):
        ref = self.ref_msm.eigenvectors_[:, 1]
        est = self.modeller.full_eigenvec

        # Only relative sign matters, pick the best.
        diff1 = ref - est
        diff2 = ref + est
        errorval = min(np.dot(diff1, diff1), np.dot(diff2, diff2))
        self.errors_over_time += [errorval]

        log.debug("Eigenvec L2 Error: %g\t Threshold: %g", errorval,
                  self.tolerance)

        return errorval < self.tolerance

    def plot(self, axs, sstate):
        top, bot = axs[0:2]

        ref = self.ref_msm.eigenvectors_[:, 1]
        est = self.modeller.full_eigenvec

        scatter_eigenvector(top, self.centers, est, ref)
        top.set_title('$||\hat{\phi}_1-\phi_1||$')

        # Bottom (over time)
        bot.plot(self.errors_over_time, 'o-')
        bot.axhline(0, c='k')
        bot.set_xlabel('Time')

    @property
    def n_plots(self):
        return 2


class TMatFro(CentroidConvergenceChecker):
    """Frobenius norm."""

    def check_convergence(self, params):
        est = self.modeller.full_tmat.todense()
        ref = self.ref_msm.transmat_.todense()

        errorval = scipy.linalg.norm(est - ref, ord='fro')
        self.errors_over_time += [errorval]

        log.debug("Frobenius Norm: %g\t Threshold: %g", errorval,
                  self.tolerance)

        return errorval < self.tolerance

    def plot(self, axs, sstate):
        top, bot = axs[0:2]

        top.set_title('Tmat Frobenius Norm')

        bot.plot(self.errors_over_time, 'o-')
        bot.axhline(0, c='k')
        bot.set_xlabel('Time')

    @property
    def n_plots(self):
        return 2

