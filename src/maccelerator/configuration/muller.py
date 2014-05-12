__author__ = 'harrigan'

from maccelerator.configuration import base_classes as bc
import maccelerator as maccel


class MullerSimulator(maccel.OpenMMSimulator):
    def __init__(self):
        # TODO: Pass in arguments
        super(maccel.OpenMMSimulator, self).__init__()

    def simulate(self, n_steps):
        super(maccel.OpenMMSimulator, self).simulate(n_steps, minimize=False,
                                                     random_initial_velocities=True)


class MullerConfiguration(bc.OpenMMConfiguration):
    def __init__(self):
        super(bc.OpenMMConfiguration, self).__init__()

        self.simulator = MullerSimulator()
