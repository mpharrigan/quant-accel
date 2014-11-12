import logging
import itertools

from mixtape.markovstatemodel import MarkovStateModel

from mixtape.cluster import KMeans
from mixtape.featurizer import DihedralFeaturizer

from ..simulate import TMatSimulator
from ..model import TMatModeller
from ..adapt import RandomAdapter
from .base import TMatConfiguration
from ..convergence import TMatConvergenceChecker
from ..param import AdaptiveParams
from ..files import get_fn


log = logging.getLogger(__name__)


def generate_alanine_msm(ala):
    """Make a small transition matrix from Alanine trajectories

    We featurize using phi / psi angles

    Note: This function is not-deterministic, although it would be useful
    if it were, so testing could be conducted.
    """

    # Featurize
    dihed = DihedralFeaturizer(['phi', 'psi'], sincos=False)
    feat_trajs = dihed.transform(ala['trajectories'])

    # Cluster
    kmeans = KMeans(n_clusters=20, random_state=52)
    kmeans.fit(feat_trajs)

    # Build MSM
    msm = MarkovStateModel(lag_time=3, verbose=False)
    msm.fit(kmeans.labels_)

    return msm, kmeans


class AlanineParams(AdaptiveParams):
    """Returned by get_param_grid which must be implemented."""

    def __init__(self, spt, tpr, adapt_lt=1, build_lt=1,
                 post_converge=10, run_id=0, step_res=8):
        super().__init__(spt=spt, tpr=tpr, adapt_lt=adapt_lt, build_lt=build_lt,
                         post_converge=post_converge, run_id=run_id,
                         step_res=step_res)


class AlanineConfiguration(TMatConfiguration):
    def __init__(self):
        super().__init__(get_fn('ala.msm.pickl'), get_fn('ala.centers.h5'))

    def defaults(config):
        # Simulator
        config.simulator_class = TMatSimulator

        # Modeller
        config.modeller_class = TMatModeller

        # Convergence checker
        config.convchecker_class = TMatConvergenceChecker

        # Adapter
        config.adapter_class = RandomAdapter

        # Define a function to yield combinations of parameters
        def get_param_grid(run_id):
            spts = [10, 100]
            tprs = [1, 10]

            for spt, tpr in itertools.product(spts, tprs):
                yield AlanineParams(spt=spt, tpr=tpr, run_id=run_id)

        config.get_param_grid = get_param_grid

        # Tolerance scaling factor
        config.tolerance_scale = 1.0



