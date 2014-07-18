"""Named hybrid convergence criteria."""

from .base import HybridConvergenceChecker
from .centroid import PopulationCentroidTVD, EigenvecCentroid, EigenvecL2, \
    TMatFro
from .general import TimescaleDistance

__author__ = 'harrigan'


class TMatConvergenceChecker(HybridConvergenceChecker):
    """A Hybrid convergence criteria for transition matrix simulations

    :param tolerance_scale: We define a default set of tolerances, which
            you can scale using this number
    """

    def __init__(self, modeller, centers, ref_msm, tolerance_scale=1):
        super().__init__(
            PopulationCentroidTVD(tolerance_scale * 0.08, modeller, centers,
                                  ref_msm),
            EigenvecCentroid(tolerance_scale * 0.01, modeller, centers,
                             ref_msm),
            EigenvecL2(tolerance_scale * 0.01, modeller, centers, ref_msm),
            TMatFro(tolerance_scale * 0.5, modeller, centers, ref_msm),
            TimescaleDistance(tolerance_scale * 20, modeller, ref_msm)
        )
