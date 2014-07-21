import pickle

import mdtraj.io
from mixtape.markovstatemodel import MarkovStateModel
from mixtape.cluster import KMeans
from mixtape.featurizer import DihedralFeaturizer
import numpy as np

from ..simulate import TMatSimulator
from ..model import TMatModeller
from ..adapt import RandomAdapter
from ..configuration import TMatConfiguration
from ..convergence.hybrid import TMatConvergenceChecker
from ..param import AdaptiveParams


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
    msm = MarkovStateModel(n_states=20, lag_time=3, n_timescales=5)
    msm.fit(kmeans.labels_)

    return msm, kmeans


class AlanineSimulator(TMatSimulator):
    pass


class AlanineModeller(TMatModeller):
    def __init__(self, tot_n_states):
        super().__init__(tot_n_states)

    def seed_state(self, params, sstate_out_fn):
        sstate = np.asarray([0] * params.tpr, dtype=int)
        mdtraj.io.saveh(sstate_out_fn, starting_states=sstate)
        return sstate

    def model(self, traj_fns, params):
        trajs = [mdtraj.io.loadh(fn, 'state_traj') for fn in traj_fns]
        return super()._model(trajs, lagtime=params.adapt_lt)


class AlanineParams(AdaptiveParams):
    @property
    def build_lt(self):
        return 1

    @property
    def adapt_lt(self):
        return 1

    @property
    def post_converge(self):
        # TODO Change
        return 10


class AlanineAdapter(RandomAdapter):
    def adapt(self, params, sstate_out_fn):
        indices = np.asarray(super()._adapt(params), dtype=int)
        mdtraj.io.saveh(sstate_out_fn, starting_states=indices)
        return indices

    @property
    def sstatefn(self):
        return 'sstate-{round_i}.h5'


class AlanineConvchecker(TMatConvergenceChecker):
    pass


class AlanineConfiguration(TMatConfiguration):
    def __init__(self, ref_msm_fn, centers_fn):
        super().__init__()

        # Load reference MSM
        with open(ref_msm_fn, 'rb') as ref_msm_f:
            ref_msm = pickle.load(ref_msm_f)

        # Load cluster centers for visualization
        centers = mdtraj.io.loadh(centers_fn, 'cluster_centers')

        # Set fields
        self.simulator = AlanineSimulator(ref_msm.transmat_)
        self.modeller = AlanineModeller(tot_n_states=ref_msm.n_states)
        self.convchecker = AlanineConvchecker(self.modeller, centers, ref_msm)
        self.adapter = AlanineAdapter(self.modeller)


