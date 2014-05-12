__author__ = 'harrigan'

import unittest
import maccelerator as maccel
from maccelerator.testing.utils import get_fn
import numpy as np
import scipy.sparse


class TestMullerSampling(unittest.TestCase):
    def setUp(self):
        self.config = maccel.configuration.muller.MullerConfiguration()

    def test_sstate(self):
        sstate = self.config.modeler.seed_state()
        np.testing.assert_array_equal(sstate.xyz, [[[0.5, 0.0, 0.0]]])

    def test_sampling_length(self):
        sstate = self.config.modeler.seed_state()

        # TODO: This has to write to disk
        traj = self.config.simulator.simulate(sstate, 100)

        self.assertEqual(traj.n_frames, 100)


class TestSimple(unittest.TestCase):
    def setUp(self):
        pass


    def test_sample(self):
        pass


class TestTMatSampling(unittest.TestCase):
    def setUp(self):
        tmat = scipy.sparse.csr_matrix(
            np.diag(np.ones(5), 1)  # 6 state linear transition model
        )

        self.simulator = maccel.TMatSimulator(tmat)

    def test_length(self):
        seq = self.simulator.simulate(sstate=0, n_steps=6)
        self.assertEqual(len(seq), 6)

    def test_states(self):
        seq = self.simulator.simulate(sstate=0, n_steps=6)
        np.testing.assert_array_equal(seq, range(6))



