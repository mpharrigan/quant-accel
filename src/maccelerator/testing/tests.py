__author__ = 'harrigan'

import unittest
import tempfile

import numpy as np
import scipy.sparse

import maccelerator as maccel


class TestMullerSampling(unittest.TestCase):
    def setUp(self):
        self.config = maccel.MullerConfiguration()

    def test_sstate(self):
        sstate = self.config.modeller.seed_state()
        for ss in sstate:
            np.testing.assert_array_equal(ss.xyz, [[[0.5, 0.0, 0.0]]])

    def test_sampling_length(self):
        sstate = self.config.modeller.seed_state()

        # TODO: This has to write to disk
        traj = self.config.simulator.simulate(sstate, 100)

        self.assertEqual(traj.n_frames, 100)


class TestSimple(unittest.TestCase):
    def setUp(self):
        self.config = maccel.SimpleConfiguration()


    def test_sstate(self):
        sstate = self.config.modeller.seed_state(8)

        self.assertEqual(len(sstate), 8)

        for ss in sstate:
            self.assertEqual(ss, 0)


    def test_sample(self):
        sstate = self.config.modeller.seed_state(1)
        seq = self.config.simulator.simulate(sstate[0], 10, traj_out_fn=None)
        np.testing.assert_array_equal(seq, range(10))


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
        configuration = maccel.SimpleConfiguration()
        param = maccel.SimpleParams(spt=10, tpr=10)
        rundir = tempfile.mkdtemp()

        self.run = maccel.MAccelRun(configuration, param, rundir)

    def test_run(self):
        self.run.run()

        self.assertEqual(len(self.run.trajs), 9)
