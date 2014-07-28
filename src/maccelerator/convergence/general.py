"""General convergence criteria."""

from .base import ConvergenceChecker, Convergence
import numpy as np

__author__ = 'harrigan'


class TimescaleDistanceConvergence(Convergence):
    def __init__(self, converged, est_timescales, ref_timescales,
                 errors_over_time):
        super().__init__(converged, errors_over_time)

        self.est_timescales = est_timescales
        self.ref_timescales = ref_timescales
        self.errors_over_time = errors_over_time

    def plot(self, axs, sstate):
        top, bot = axs[0:2]

        est = self.est_timescales
        est = est[~np.isnan(est)]
        est = est[est > 0]
        ref = self.ref_timescales

        top.set_title('Timescale difference')
        top.hlines(est, 0, 0.5, 'r', label='Est')
        top.hlines(ref, 0.5, 1, 'b', label='Ref')
        top.legend(loc='best')
        top.set_yscale('log')

        self._plot_bottom(bot)

    @property
    def n_plots(self):
        return 2


class TimescaleDistance(ConvergenceChecker):
    def __init__(self, tolerance, ref_msm, n_timescales=1):
        super().__init__(tolerance)
        self.ref_msm = ref_msm

    def check_convergence(self, model, params):
        est_timescales = model.timescales / model.lagtime
        ref_timescales = self.ref_msm.timescales_ / self.ref_msm.lag_time

        # TODO: Take the difference of the logs

        # TODO: Include more than one timescale
        est = est_timescales[0]
        ref = ref_timescales[0]

        errorval = np.abs(est - ref)
        self.errors_over_time += [errorval]
        converged = errorval < self.tolerance
        return TimescaleDistanceConvergence(
            converged, est_timescales, ref_timescales, self.errors_over_time
        )


