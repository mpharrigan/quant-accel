import pickle

import mdtraj.io

from ..simulate import TMatSimulator
from ..model import TMatModeller
from ..adapt import RandomAdapter
from ..configuration import TMatConfiguration
from ..convergence.hybrid import TMatConvergenceChecker
from ..param import AdaptiveParams


# TODO: Put make reference function in here

class AlanineSimulator(TMatSimulator):
    pass


class AlanineModeller(TMatModeller):
    def __init__(self, tot_n_states):
        super().__init__(tot_n_states)

    def seed_state(self, params):
        return [0] * params.tpr

    def model(self, traj_fns, params):
        trajs = [mdtraj.io.loadh(fn, 'state_traj') for fn in traj_fns]
        super()._model(trajs, lagtime=params.adapt_lt)


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
        return 1

    @property
    def threshold(self):
        return 0.05


class AlanineAdapter(RandomAdapter):
    pass


class AlanineConvchecker(TMatConvergenceChecker):
    pass


class AlanineConfiguration(TMatConfiguration):
    def __init__(self, ref_msm_fn, centers_fn):
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


