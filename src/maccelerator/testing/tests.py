__author__ = 'harrigan'

import unittest

import maccelerator as maccel
import numpy as np
import scipy.sparse


class TestMullerSampling(unittest.TestCase):
    def setUp(self):
        self.config = maccel.MullerConfiguration()

    def test_sstate(self):
        sstate = self.config.modeller.seed_state()
        np.testing.assert_array_equal(sstate.xyz, [[[0.5, 0.0, 0.0]]])

    def test_sampling_length(self):
        sstate = self.config.modeller.seed_state()

        # TODO: This has to write to disk
        traj = self.config.simulator.simulate(sstate, 100)

        self.assertEqual(traj.n_frames, 100)


class TestSimple(unittest.TestCase):
    def setUp(self):
        self.config = maccel.SimpleConfiguration()


    def test_sstate(self):
        sstate = self.config.modeller.seed_state()
        self.assertEqual(sstate, 0)


    def test_sample(self):
        sstate = self.config.modeller.seed_state()
        seq = self.config.simulator.simulate(sstate, 10)


class TestTMatSampling(unittest.TestCase):
    def setUp(self):
        tmat = scipy.sparse.csr_matrix(
            np.diag(np.ones(5), 1)  # 6 state linear transition model
        )

        # TODO: Don't test like this
        self.simulator = maccel.TMatSimulator(tmat)

    def test_length(self):
        seq = self.simulator.simulate(sstate=0, n_steps=6)
        self.assertEqual(len(seq), 6)

    def test_states(self):
        seq = self.simulator.simulate(sstate=0, n_steps=6)
        np.testing.assert_array_equal(seq, range(6))


class TestRun(unittest.TestCase):
    def setUp(self):
        self.configuration = maccel.SimpleConfiguration()
        self.param = maccel.SimpleParams(spt=10, tpr=10)
        self.run = maccel.MAccelRun(self.configuration, self.param)

    def test_run(self):
        self.run.run()
