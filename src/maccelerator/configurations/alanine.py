from ..simulate import TMatSimulator
from ..model import TMatModeller, SortCountsAdapter
from ..configuration import TMatConfiguration
from ..check_convergence import PopulationCentroidTVD, EigenvecCentroid, \
    HybridConvergenceChecker
from ..param import AdaptiveParams

import pickle
import mdtraj.io


class AlanineSimulator(TMatSimulator):
    pass


class AlanineModeller(TMatModeller):
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
        return 0


class AlanineAdapter(SortCountsAdapter):
    pass


class AlanineConvchecker(HybridConvergenceChecker):
    def __init__(self, modeller, centers):
        super().__init__(PopulationCentroidTVD(modeller, centers),
                         EigenvecCentroid(modeller, centers))
        self.do_plots = True


class AlanineConfiguration(TMatConfiguration):
    def __init__(self, ref_msm_fn, centers_fn):
        with open(ref_msm_fn, 'rb') as ref_msm_f:
            ref_msm = pickle.load(ref_msm_f)

        centers = mdtraj.io.loadh(centers_fn, 'cluster_centers')

        super().__init__(ref_msm)
        self.simulator = AlanineSimulator(self.tmat)
        self.modeller = AlanineModeller()
        self.convchecker = AlanineConvchecker(self.modeller, centers)

        self.adapter = AlanineAdapter(self.modeller)


