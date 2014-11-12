"""Test sampling from alanine transition matrix."""

__author__ = 'harrigan'

from unittest import TestCase
from os.path import join as pjoin
import os
import unittest
import logging

import mdtraj.io
import maccelerator as maccel
from maccelerator.adapt import SStates
from maccelerator.model import Model
from maccelerator.param import AdaptiveParams
from maccelerator.testing.test_utils import get_folder
from mixtape.datasets.alanine_dipeptide import fetch_alanine_dipeptide
from maccelerator.configurations.alanine import generate_alanine_msm


logging.captureWarnings(True)
logging.disable(logging.WARNING)


class TestAlaninePrep(TestCase):
    def test_make(self):
        ala = fetch_alanine_dipeptide()
        msm, kmeans = generate_alanine_msm(ala)

        self.assertEqual(len(kmeans.labels_), 20)
        self.assertEqual(msm.transmat_.shape[1], 20)


class TestRun(TestCase):
    """Test features related to doing one adaptive run."""

    do_parallel = True

    def setUp(self):
        configuration = maccel.AlanineConfiguration().apply_configuration()
        self.tpr = 2
        self.spt = 1000
        param = maccel.AlanineParams(spt=self.spt, tpr=self.tpr)
        self.rundir = get_folder('ala')
        self.run = maccel.MAccelRun(configuration, param, self.rundir,
                                    parallel=self.do_parallel)
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

    def test_partial_load(self):
        trajs = self.run.config.modeller.load_trajs(self.run.le_than(0), 555)
        self.assertEqual(len(trajs), 3)
        for i, traj in enumerate(trajs):
            with self.subTest(traj_i=i):
                self.assertEqual(len(traj), 555)

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


class TestRunNoParallel(TestRun):
    def setUp(self):
        self.do_parallel = False
        super().setUp()


if __name__ == "__main__":
    unittest.main()
