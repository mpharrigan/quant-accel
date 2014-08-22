import pickle
import logging

import mdtraj.io
from mixtape.markovstatemodel import MarkovStateModel
from mixtape.cluster import KMeans
from mixtape.featurizer import DihedralFeaturizer

from ..simulate import TMatSimulator
from ..model import TMatModeller
from ..adapt import RandomAdapter, SStates
from .base import TMatConfiguration
from ..convergence.hybrid import TMatConvergenceChecker
from ..param import AdaptiveParams

log = logging.getLogger(__name__)

_AlANINE_TEMPLATE = """
class MyAlanineConfiguration(maccel.AlanineConfiguration):
    def __init__(self):
        super().__init__(maccel.get_fn('ala.msm.pickl'),
                         maccel.get_fn('ala.centers.h5'))

    def get_param_grid(self, run_id):
        spts = [10, 100]
        tprs = [10, 100]

        for spt, tpr in itertools.product(spts, tprs):
            yield maccel.AlanineParams(spt=spt, tpr=tpr, run_id=run_id)
"""


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


class AlanineSimulator(TMatSimulator):
    pass


class AlanineModeller(TMatModeller):
    def load_trajs(self, traj_fns, up_to=None):
        if up_to is not None:
            trajs = [mdtraj.io.loadh(fn, 'state_traj')[:up_to] for fn in
                     traj_fns]
        else:
            trajs = [mdtraj.io.loadh(fn, 'state_traj') for fn in traj_fns]
        return trajs

    def lagtime(self, params):
        return params.adapt_lt


class AlanineAdapter(RandomAdapter):
    def seed_states(self, params):
        indices = [0] * params.tpr
        return SStates(indices)


class AlanineConvchecker(TMatConvergenceChecker):
    pass


class AlanineParams(AdaptiveParams):
    @property
    def build_lt(self):
        return 1

    @property
    def adapt_lt(self):
        return 1

    @property
    def post_converge(self):
        return 10


class AlanineConfiguration(TMatConfiguration):
    @classmethod
    def get_template(cls, grid_manager_name):
        fmt_dict = dict(grid_manager=grid_manager_name,
                        other_config=_AlANINE_TEMPLATE)
        temp = super().get_template(grid_manager_name)
        return temp.format(**fmt_dict)

    def __init__(self, ref_msm_fn, centers_fn):
        super().__init__()

        # Load reference MSM
        with open(ref_msm_fn, 'rb') as ref_msm_f:
            ref_msm = pickle.load(ref_msm_f)

        # Load cluster centers for visualization
        centers = mdtraj.io.loadh(centers_fn, 'cluster_centers')

        # Set fields
        self.simulator = AlanineSimulator(ref_msm.transmat_)
        self.modeller = AlanineModeller(tot_n_states=ref_msm.n_states_)
        self.convchecker = AlanineConvchecker(centers, ref_msm)
        self.adapter = AlanineAdapter()

    def get_param_grid(self, run_id=0):
        raise NotImplementedError()


