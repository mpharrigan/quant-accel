"""Test sampling from alanine transition matrix."""

__author__ = 'harrigan'

from unittest import TestCase
from os.path import join as pjoin
import os
import unittest
import logging
import pickle

import mdtraj.io
import scipy.io

import maccelerator as maccel
from maccelerator.adapt import SStates
from maccelerator.model import Model
from maccelerator.param import AdaptiveParams
from maccelerator.testing.utils import get_folder, get_fn
from maccelerator.testing.make_reference_data import make_alanine_reference_data


# Disable logging during test
logging.captureWarnings(True)
logging.disable(logging.WARNING)


class TestAlaninePrep(TestCase):
    def setUp(self):
        self.rundir = get_folder('ala_prep')

    @unittest.skip
    def test_make(self):
        dirname = self.rundir
        make_alanine_reference_data(dirname)
        # Save cluster centers
        centers = mdtraj.io.loadh(pjoin(dirname, 'ala.centers.h5'),
                                  'cluster_centers')

        self.assertEqual(len(centers), 20)

        # Save transition matrix
        tmat = scipy.io.mmread(pjoin(dirname, 'ala.mtx'))

        self.assertEqual(tmat.shape[1], 20)

        # Save MSM Object
        with open(pjoin(dirname, 'ala.msm.pickl'), 'rb') as f:
            msm = pickle.load(f)

        self.assertIsNotNone(msm)


class TestRun(TestCase):
    """Test features related to doing one adaptive run."""

    def setUp(self):
        configuration = maccel.AlanineConfiguration(get_fn('ala.msm.pickl'),
                                                    get_fn('ala.centers.h5'))
        self.tpr = 3
        self.spt = 1000
        param = maccel.AlanineParams(spt=self.spt, tpr=self.tpr)
        self.rundir = get_folder('ala')
        self.run = maccel.MAccelRun(configuration, param, self.rundir)
        self.run.run()

    def test_num_trajs(self):
        # 2 trajectories per round
        self.assertEqual(len(self.run.trajs), self.run.n_rounds)

        for round_i, trajs in self.run.trajs.items():
            with self.subTest(round_i=round_i):
                self.assertEqual(len(trajs), self.tpr)

    def test_param_file(self):
        paramfn = pjoin(self.rundir, 'params.pickl')
        params = AdaptiveParams.load(paramfn)

        self.assertEqual(self.tpr, params.tpr)
        self.assertEqual(self.spt, params.spt)

    def test_run_file(self):
        runfn = pjoin(self.rundir, 'run.pickl')
        run = maccel.MAccelRun.load(runfn)

        figfn = 'plot-{:04d}'.format(1)

        self.assertEqual(run.rundir, self.rundir)
        self.assertEqual(run.params.tpr, self.tpr)
        self.assertEqual(run.config.file.plot_fn(1),
                         pjoin(self.rundir, 'figs', figfn))


    def test_trajectory_files(self):
        # Make sure files exist where they should and don't where they shouldn't

        nrounds = len(self.run.trajs)

        for round_i in range(nrounds * 2):
            for traj_i in range(self.tpr * 2):
                should_exist = round_i < nrounds and traj_i < self.tpr
                with self.subTest(round_i=round_i, traj_i=traj_i,
                                  should_exist=should_exist):
                    trajoutfn = 'round-{round_i}/traj-{traj_i}.h5'.format(
                        round_i=round_i, traj_i=traj_i)
                    trajoutfn = pjoin(self.rundir, 'trajs', trajoutfn)
                    if should_exist:
                        self.assertTrue(os.path.exists(trajoutfn))
                        traj = mdtraj.io.loadh(trajoutfn, 'state_traj')
                        self.assertEqual(len(traj), self.spt)
                    else:
                        self.assertFalse(os.path.exists(trajoutfn))

    def test_sstate_files(self):
        nrounds = len(self.run.trajs)
        for round_i in range(nrounds * 2):
            should_exist = round_i <= nrounds
            with self.subTest(round_i=round_i, should_exist=should_exist):
                sstatefn = 'sstate-{round_i}.pickl'.format(round_i=round_i)
                sstatefn = pjoin(self.rundir, 'sstates', sstatefn)
                if should_exist:
                    self.assertTrue(os.path.exists(sstatefn))
                    sstate = SStates.load(sstatefn)
                    self.assertEqual(len(sstate.items()), self.tpr)
                else:
                    self.assertFalse(os.path.exists(sstatefn))

    def test_msm_files(self):
        nrounds = len(self.run.trajs)
        for round_i in range(nrounds * 2):
            should_exist = round_i < nrounds
            with self.subTest(round_i=round_i, should_exist=should_exist):
                msmfn = 'msm-{round_i}.pickl'.format(round_i=round_i)
                msmfn = pjoin(self.rundir, 'msms', msmfn)
                if should_exist:
                    self.assertTrue(os.path.exists(msmfn))
                    model = Model.load(msmfn)
                    self.assertEqual(model.tot_n_states, 20)
                else:
                    self.assertFalse(os.path.exists(msmfn))


if __name__ == "__main__":
    unittest.main()
