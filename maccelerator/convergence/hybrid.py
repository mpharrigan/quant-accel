"""Named hybrid convergence criteria."""

from .base import SupConvergenceChecker
from .centroid import PopulationCentroidTVD, EigenvecCentroid, EigenvecL2, \
    TMatFro
from .general import TimescaleDistance

__author__ = 'harrigan'


class TMatConvergenceChecker(SupConvergenceChecker):
    """A Hybrid convergence criteria for transition matrix simulations
    """

    @classmethod
    def get_sub_checkers(cls, config):
        centers = config.centers
        ref_msm = config.ref_msm
        tolerance_scale = config.tolerance_scale
        return [
            PopulationCentroidTVD(tolerance_scale * 0.08, centers, ref_msm),
            EigenvecCentroid(tolerance_scale * 0.01, centers, ref_msm),
            EigenvecL2(tolerance_scale * 0.08, centers, ref_msm),
            TMatFro(tolerance_scale * 0.5, centers, ref_msm),
            TimescaleDistance(tolerance_scale * 0.2, ref_msm)
        ]
