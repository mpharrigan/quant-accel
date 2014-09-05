import logging
import itertools

from mixtape.markovstatemodel import MarkovStateModel

from mixtape.cluster import KMeans
from mixtape.featurizer import DihedralFeaturizer

from ..simulate import TMatSimulator

from ..model import TMatModeller
from ..adapt import RandomAdapter
from .base import TMatConfiguration
from ..convergence.hybrid import TMatConvergenceChecker
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
    msm = MarkovStateModel(lag_time=3)
    msm.fit(kmeans.labels_)

    return msm, kmeans


class AlanineParams(AdaptiveParams):
    """Returned by get_param_grid which must be implemented."""
    pass


class AlanineConfiguration(TMatConfiguration):
    def __init__(self):
        super().__init__(get_fn('ala.msm.pickl'), get_fn('ala.centers.h5'))

    def defaults(config):
        # Lag time to build the model
        setattr(AlanineParams, 'build_lt', property(lambda self: 1))

        # Lag time to generate counts for adapting
        setattr(AlanineParams, 'adapt_lt', property(lambda self: 1))

        # Number of rounds to do after convergence is reached
        setattr(AlanineParams, 'post_converge', property(lambda self: 10))

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


    def apply_configuration(self):
        # TODO: Make all of these things just take in "config" and duck

        self.simulator = self.simulator_class(self.ref_msm.transmat_)
        self.modeller = self.modeller_class(tot_n_states=self.ref_msm.n_states_)
        self.convchecker = self.convchecker_class(self.centers, self.ref_msm)
        self.adapter = self.adapter_class()


