"""Test sampling from alanine transition matrix."""

__author__ = 'harrigan'

from unittest import TestCase
from os.path import join as pjoin
import os
import unittest
import logging

import numpy as np
from numpy.testing import assert_array_equal
import mdtraj.io

import maccelerator as maccel
from maccelerator.testing.utils import get_folder, get_fn



# Disable logging during test
logging.captureWarnings(True)
logging.disable(logging.WARNING)


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
        for round_i, trajs in self.run.trajs.items():
            with self.subTest(round_i=round_i):
                self.assertEqual(len(trajs), self.tpr)

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
                sstatefn = 'sstate-{round_i}.h5'.format(round_i=round_i)
                sstatefn = pjoin(self.rundir, 'sstates', sstatefn)
                if should_exist:
                    self.assertTrue(os.path.exists(sstatefn))
                    sstate = mdtraj.io.loadh(sstatefn, 'starting_states')
                    self.assertEqual(len(sstate), self.tpr)
                else:
                    self.assertFalse(os.path.exists(sstatefn))

    @unittest.skip
    def test_msm_files(self):
        return False
