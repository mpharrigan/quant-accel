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

    def simulate(self, sstate, n_steps):
        return range(sstate, n_steps)


class SimpleModeller(Modeller):
    def __init__(self):
        super().__init__()

    def seed_state(self):
        return 0


class SimpleAdapter(Adapter):
    def __init__(self, modeller):
        super().__init__()
        self.modeller = modeller

    def adapt(self, n_tpr):
        counts = self.modeller.counts
        return np.array([np.argmin(counts)] * n_tpr)


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
