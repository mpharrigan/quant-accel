"""Named hybrid convergence criteria."""

from .base import HybridConvergenceChecker
from .centroid import PopulationCentroidTVD, EigenvecCentroid, EigenvecL2, \
    TMatFro
from .general import TimescaleDistance

__author__ = 'harrigan'


class TMatConvergenceChecker(HybridConvergenceChecker):
    def __init__(self, modeller, centers, ref_msm):
        super().__init__(
            PopulationCentroidTVD(modeller, centers, ref_msm),
            EigenvecCentroid(modeller, centers, ref_msm),
            EigenvecL2(modeller, centers, ref_msm),
            TMatFro(modeller, centers, ref_msm),
            TimescaleDistance(modeller, ref_msm)

        )
