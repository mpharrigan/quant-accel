"""Tests for the overarching framework of grid.py and run.py

Also test the 'simple' configuration that is included for testing.
"""

__author__ = 'harrigan'

import os
from os.path import join as pjoin
import unittest
from unittest import TestCase
import logging

import numpy as np
from numpy.testing import assert_array_equal

from maccelerator.testing.test_utils import get_folder
import maccelerator as maccel



# Disable logging during test
logging.disable(logging.WARNING)


class TestSimpleStartingStates(TestCase):
    """Make sure starting state generation is working."""

    def setUp(self):
        self.config = maccel.SimpleConfiguration().apply_configuration()
        params = maccel.configurations.SimpleParams(spt=10, tpr=8)
        rundir = get_folder('s1')
        self.out_fn = pjoin(rundir, 'test_sstate')
        self.sstate = self.config.seed_state(params)
        self.sstate.save(self.out_fn)

    def test_length(self):
        self.assertEqual(len(self.sstate.items()), 8)

    def test_value(self):
        for ss in self.sstate.items():
            with self.subTest(ss=ss):
                self.assertEqual(ss, 0)

    def test_save_sstates(self):
        full_ofn = '{}.pickl'.format(self.out_fn)
        self.assertTrue(os.path.exists(full_ofn))

        # TODO: Load


class TestSimpleSample(TestCase):
    """Make sure sampling (trajectory generation) is working."""

    def setUp(self):
        self.config = maccel.SimpleConfiguration().apply_configuration()
        params = maccel.configurations.SimpleParams(spt=10, tpr=2)
        rundir = get_folder('s2')
        self.out_fn = pjoin(rundir, 'test_sample.npy')
        self.traj = self.config.simulate(0, params, traj_out_fn=self.out_fn)

    def test_save_trajectories(self):
        self.assertTrue(os.path.exists(self.out_fn))

        traj = np.load(self.out_fn)
        assert_array_equal(traj, np.arange(10))


class TestRun(TestCase):
    """Test features related to doing one adaptive run."""

    do_parallel = False

    def setUp(self):
        configuration = maccel.SimpleConfiguration().apply_configuration()
        param = maccel.SimpleParams(spt=10, tpr=10)
        rundir = get_folder('s3')
        self.rundir = rundir

        self.run = maccel.MAccelRun(configuration, param, rundir,
                                    parallel=self.do_parallel)
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
                        traj = np.load(trajoutfn)
                        assert_array_equal(traj, 9 * round_i + np.arange(10))
                    else:
                        self.assertFalse(os.path.exists(trajoutfn))

    def test_sstate_files(self):
        for round_i in range(8 * 2):
            should_exist = round_i <= 8
            with self.subTest(round_i=round_i, should_exist=should_exist):
                sstatefn = 'sstate-{round_i}.pickl'.format(round_i=round_i)
                sstatefn = pjoin(self.rundir, 'sstates', sstatefn)
                if should_exist:
                    self.assertTrue(os.path.exists(sstatefn))
                    sstate = np.load(sstatefn)
                    assert_array_equal(sstate.items(),
                                       9 * round_i + np.zeros(10))
                else:
                    self.assertFalse(os.path.exists(sstatefn))


class TestRunNoParallel(TestRun):
    def setUp(self):
        self.do_parallel = False
        super().setUp()


class TestGrid(TestCase):
    """Test features related to doing a collection of adaptive runs."""

    do_parallel = True

    def setUp(self):
        configuration = maccel.SimpleConfiguration().apply_configuration()
        self.griddir = get_folder('s4')
        self.grid = maccel.MAccelGrid(configuration, self.griddir,
                                      parallel=self.do_parallel)
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
                self.assertEqual(len(os.listdir(pjoin(self.griddir, od))), 7)

    def test_griddir_attribute(self):
        self.assertEqual(self.griddir, self.grid.griddir)

    def test_num_folders(self):
        self.assertEqual(len(os.listdir(self.griddir)), 2)

    def test_iterable_paramlist(self):
        for p in self.grid.config.get_param_grid():
            with self.subTest():
                self.assertTrue(hasattr(p, 'spt'))


class TestGridNoParallel(TestGrid):
    def setUp(self):
        self.do_parallel = False
        super().setUp()


if __name__ == "__main__":
    unittest.main()