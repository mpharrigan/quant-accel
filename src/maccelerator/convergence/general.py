"""General convergence criteria."""

from .base import ConvergenceChecker
import numpy as np

__author__ = 'harrigan'


class TimescaleDistance(ConvergenceChecker):
    def __init__(self, modeller, ref_msm):
        super().__init__()
        self.modeller = modeller
        self.ref_msm = ref_msm

    def check_convergence(self, params):
        est = self.modeller.msm.timescales_[0] / self.modeller.msm.lag_time
        ref = self.ref_msm.timescales_[0] / self.ref_msm.lag_time

        errorval = est - ref
        self.errors_over_time += [errorval]
        return np.abs(errorval) < params.threshold

    def plot(self, axs, sstate):
        top, bot = axs[0:2]

        top.set_title('Timescale difference')

        bot.plot(self.errors_over_time, 'o-')
        bot.axhline(0, c='k')
        bot.set_xlabel('Time')

    @property
    def n_plots(self):
        return 2

