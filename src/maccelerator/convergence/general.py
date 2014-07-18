"""General convergence criteria."""

from .base import ConvergenceChecker
import numpy as np

__author__ = 'harrigan'


class TimescaleDistance(ConvergenceChecker):
    def __init__(self, tolerance, modeller, ref_msm):
        super().__init__(tolerance)
        self.modeller = modeller
        self.ref_msm = ref_msm

    def check_convergence(self, params):
        est = self.modeller.msm.timescales_[0] / self.modeller.msm.lag_time
        ref = self.ref_msm.timescales_[0] / self.ref_msm.lag_time

        errorval = np.abs(est - ref)
        self.errors_over_time += [errorval]
        return errorval < self.tolerance

    def plot(self, axs, sstate):
        top, bot = axs[0:2]

        est = self.modeller.msm.timescales_ / self.modeller.msm.lag_time
        est = est[~np.isnan(est)]
        est = est[est > 0]
        ref = self.ref_msm.timescales_ / self.ref_msm.lag_time

        top.set_title('Timescale difference')
        top.hlines(est, 0, 0.5, 'r', label='Est')
        top.hlines(ref, 0.5, 1, 'b', label='Ref')
        top.legend(loc='best')
        top.set_yscale('log')

        bot.plot(self.errors_over_time, 'o-')
        bot.axhline(0, c='k')
        bot.set_xlabel('Time')

    @property
    def n_plots(self):
        return 2

