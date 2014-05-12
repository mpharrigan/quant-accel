"""Test openmm."""

import unittest
import maccelerator as maccel
from maccelerator.testing.utils import get_fn

#TODO: Test configurations *as well*
class TestOpenMMSampling(unittest.TestCase):

    def setUp(self):
        self.simulator = maccel.OpenMMSimulator(report_stride=10)
        self.simulator.deserialize(get_fn('muller_sys.xml'), get_fn('muller_int.xml'))

    def test_sample(self):
        seq = self.simulator.simulate(sstate, 100)


if __name__ == "__main__":
    unittest.main()
