"""A configuration that doesn't actually do anything,
but can be used for testing.
"""

from . import base_classes as bc
from ..simulate import Simulator
from ..model import Modeller


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


class SimpleConfiguration(bc.Configuration):
    def __init__(self):
        super().__init__()

        self.modeller = SimpleModeller()
        self.simulator = SimpleSimulator()
