"""Convergence criteria that depend on having a consistent state definition."""

import logging as log

import numpy as np
from matplotlib import pyplot as plt
import scipy.sparse
import scipy.linalg

from .base import ConvergenceChecker
from .base import distribution_norm_tvd


__author__ = 'harrigan'


class CentroidConvergenceChecker(ConvergenceChecker):
    """For when we have consistent state decompositions."""

    def __init__(self, modeller, centers, ref_msm):
        super().__init__()

        self.modeller = modeller
        self.centers = centers
        self.ref_msm = ref_msm


class PopulationCentroidTVD(CentroidConvergenceChecker):
    """Compare equilibrium populations of states to a reference MSM."""

    def __init__(self, modeller, centers, ref_msm):
        super().__init__(modeller, centers, ref_msm)
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
        bot.axhline(0, c='k')
        bot.set_xlabel('Time')

    @property
    def n_plots(self):
        return 2


class EigenvecCentroid(CentroidConvergenceChecker):
    """Compare the first eigenvector to that of a reference MSM."""

    def check_convergence(self, params):
        ref = self.ref_msm.eigenvectors_[:, 1]
        refs = scipy.sparse.csr_matrix(ref)
        tmat = self.modeller.full_tmat

        ref_val = self.ref_msm.eigenvalues_[1]

        # <y, Poy> / <y, y>
        est_val = refs.dot(refs.dot(tmat).T)[0, 0] / np.dot(ref, ref)

        errorval = ref_val - est_val
        self.errors_over_time += [errorval]

        log.debug("Eigenvec Error: %g\t Threshold: %g", errorval,
                  params.threshold)

        return errorval < params.threshold


    def plot(self, axs, sstate):
        top, bot = axs[0:2]

        est = self.modeller.full_eigenvec
        ref = self.ref_msm.eigenvectors_[:, 1]
        scale = 2e3

        rpos = np.where(ref >= 0)[0]
        rneg = np.where(ref < 0)[0]
        epos = np.where(est >= 0)[0]
        eneg = np.where(est < 0)[0]

        # Reference
        top.scatter(self.centers[rpos, 0], self.centers[rpos, 1], linewidths=0,
                    s=scale * ref[rpos], c='coral')
        top.scatter(self.centers[rneg, 0], self.centers[rneg, 1], linewidths=0,
                    s=-scale * ref[rneg], c='skyblue')

        # Estimate
        top.scatter(self.centers[epos, 0], self.centers[epos, 1], linewidths=2,
                    facecolors='none', edgecolors='r', s=scale * est[epos])
        top.scatter(self.centers[eneg, 0], self.centers[eneg, 1], linewidths=2,
                    facecolors='none', edgecolors='b', s=-scale * est[eneg])

        top.set_title('1st Eigenvec')


        # Bottom (over time)
        bot.plot(self.errors_over_time, 'o-')
        bot.axhline(0, c='k')
        bot.set_xlabel('Time')

    @property
    def n_plots(self):
        return 2


class EigenvecL2(EigenvecCentroid):
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
                  params.threshold)

        return errorval < params.threshold


class TMatFro(CentroidConvergenceChecker):
    """Frobenius norm."""

    def check_convergence(self, params):
        est = self.modeller.full_tmat.todense()
        ref = self.ref_msm.transmat_.todense()

        errorval = scipy.linalg.norm(est - ref, ord='fro')
        self.errors_over_time += [errorval]

        log.debug("Frobenius Norm: %g\t Threshold: %g", errorval,
                  params.threshold)

        return errorval < params.threshold

    def plot(self, axs, sstate):
        top, bot = axs[0:2]

        bot.plot(self.errors_over_time, 'o-')
        bot.axhline(0, c='k')

    @property
    def n_plots(self):
        return 2

