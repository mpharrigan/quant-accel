"""A configuration that doesn't actually do anything,
but can be used for testing.
"""

from ..simulate import Simulator
from ..model import Modeller
from ..configuration import Configuration


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


class SimpleConfiguration(Configuration):
    def __init__(self):
        super().__init__()

        self.modeller = SimpleModeller()
        self.simulator = SimpleSimulator()
