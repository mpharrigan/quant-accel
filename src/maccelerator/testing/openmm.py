"""Test openmm."""

import unittest
import maccelerator as maccel

class TestOpenMMSampling(unittest.TestCase):

    def setUp(self):
        pass

    def test_sample(self):
        seq = self.openmm.sample(5)
