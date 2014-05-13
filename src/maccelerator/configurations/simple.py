"""A configuration that doesn't actually do anything,
but can be used for testing.
"""

import numpy as np

from ..simulate import Simulator
from ..model import Modeller, Adapter
from ..configuration import Configuration
from ..param import AdaptiveParams


class SimpleSimulator(Simulator):
    def __init__(self):
        super().__init__()

    def simulate(self, sstate, n_steps, traj_out_fn):
        traj = range(sstate, n_steps)

        if traj_out_fn is not None:
            np.save(traj_out_fn, traj)
        return traj

    @property
    def trajfn(self):
        return "traj-{traj_i}.npy"


class SimpleModeller(Modeller):
    def __init__(self):
        super().__init__()
        self.trajs = []

    def model(self, traj_fns):
        self.trajs = [np.load(tfn) for tfn in traj_fns]

    def check_convergence(self):
        return len(self.trajs) > 40

    def seed_state(self, tpr):
        return [0] * tpr


class SimpleAdapter(Adapter):
    def __init__(self, modeller):
        super().__init__()
        self.modeller = modeller

    def adapt(self, n_tpr):
        ret = np.array([np.max(t) for t in self.modeller.trajs[-n_tpr:]])
        assert len(ret) == n_tpr
        return ret


class SimpleParams(AdaptiveParams):
    def __init__(self, spt, tpr):
        super().__init__(spt, tpr, run_id=0)

    @property
    def post_converge(self):
        return 5


class SimpleConfiguration(Configuration):
    def __init__(self):
        super().__init__()

        self.modeller = SimpleModeller()
        self.adapter = SimpleAdapter(self.modeller)
        self.simulator = SimpleSimulator()
