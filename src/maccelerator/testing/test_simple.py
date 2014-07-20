"""Tests for the overarching framework of grid.py and run.py

Also test the 'simple' configuration that is included for testing.
"""

__author__ = 'harrigan'

import tempfile
import os
from os.path import join as pjoin
from unittest import TestCase
import logging

import numpy as np

import maccelerator as maccel


# Disable logging during test
logging.disable(logging.WARNING)


class TestSimpleStartingStates(TestCase):
    """Make sure starting state generation is working."""

    def setUp(self):
        self.config = maccel.SimpleConfiguration()
        params = maccel.configurations.SimpleParams(spt=10, tpr=8)
        self.sstate = self.config.modeller.seed_state(params)

    def test_length(self):
        self.assertEqual(len(self.sstate), 8)

    def test_value(self):
        for ss in self.sstate:
            with self.subTest(ss=ss):
                self.assertEqual(ss, 0)

    def test_save(self):
        # TODO
        self.assertTrue(True)


class TestSimpleSample(TestCase):
    """Make sure sampling (trajectory generation) is working."""

    def setUp(self):
        self.config = maccel.SimpleConfiguration()
        params = maccel.configurations.SimpleParams(spt=10, tpr=2)

        self.seqs = self.config.simulator.simulate(0, params,
                                                   traj_out_fn=None)

    def test_num_trajectories(self):
        self.assertEqual(len(self.seqs), 2)

    def test_value_trajectories(self):
        for i, seq in enumerate(self.seqs):
            with self.subTest(i=i):
                self.assertTrue(np.array_equal(seq, np.arange(10)))

    def test_save_trajectories(self):
        # TODO
        self.assertTrue(True)


class TestRun(TestCase):
    """Test features related to doing one adaptive run."""

    def setUp(self):
        configuration = maccel.SimpleConfiguration()
        param = maccel.SimpleParams(spt=10, tpr=10)
        rundir = tempfile.mkdtemp()
        self.rundir = rundir

        self.run = maccel.MAccelRun(configuration, param, rundir)
        self.run.run()

    def test_num_rounds(self):
        # Converges in 4, post converge is 4
        self.assertEqual(len(self.run.trajs), 4 + 4)

    def test_num_trajs(self):
        # 10 trajectories per round
        for round_i, trajs in self.run.trajs.items():
            with self.subTest(round_i=round_i):
                self.assertEqual(len(trajs), 10)

    def test_trajectory_files(self):
        # Make sure files exist where they should and don't where they shouldn't
        for round_i in range(8 * 2):
            for traj_i in range(10 * 2):
                should_exist = round_i < 8 and traj_i < 10
                with self.subTest(round_i=round_i, traj_i=traj_i,
                                  should_exist=should_exist):
                    trajoutfn = 'round-{round_i}/traj-{traj_i}.npy'.format(
                        round_i=round_i, traj_i=traj_i)
                    trajoutfn = pjoin(self.rundir, 'trajs', trajoutfn)
                    if should_exist:
                        self.assertTrue(os.path.exists(trajoutfn))
                    else:
                        self.assertFalse(os.path.exists(trajoutfn))


class TestGrid(TestCase):
    """Test features related to doing a collection of adaptive runs."""

    def setUp(self):
        configuration = maccel.SimpleConfiguration()
        self.griddir = tempfile.mkdtemp()
        self.grid = maccel.MAccelGrid(configuration, self.griddir)
        self.grid.grid()

    def test_name_folders(self):
        outdirs = [os.path.basename(d) for d in os.listdir(self.grid.griddir)]
        should_be = ['blt-1_alt-1_spt-10_tpr-10', 'blt-1_alt-1_spt-20_tpr-10']

        for sb in should_be:
            with self.subTest():
                self.assertTrue(sb in outdirs)

    def test_num_subfolders(self):
        outdirs = os.listdir(self.griddir)
        for od in outdirs:
            with self.subTest(outdir=od):
                # See if it has the right number of subfolders
                # Figs/, msms/, sstates/, trajs/
                self.assertTrue(len(os.listdir(pjoin(self.griddir, od))), 4)

    def test_griddir_attribute(self):
        self.assertEqual(self.griddir, self.grid.griddir)

    def test_num_folders(self):
        self.assertEqual(len(os.listdir(self.griddir)), 2)

    def test_iterable_paramlist(self):
        for p in self.grid.config.get_param_grid():
            with self.subTest():
                self.assertTrue(hasattr(p, 'spt'))
