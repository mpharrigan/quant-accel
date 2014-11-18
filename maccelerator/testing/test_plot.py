__author__ = 'harrigan'
from unittest import TestCase
from os.path import join as pjoin
import os
import logging
import unittest
import shutil

import maccelerator as maccel

from maccelerator.testing.utils import get_folder


logging.captureWarnings(True)
logging.disable(logging.WARNING)


class TestPlot(TestCase):
    do_parallel = True

    def setUp(self):
        configuration = maccel.AlanineConfiguration().apply_configuration()
        self.tpr = 3
        self.spt = 4000
        param = maccel.AlanineParams(spt=self.spt, tpr=self.tpr, step_res=-1,
                                     post_converge=2)
        self.rundir = get_folder('plot')
        run = maccel.MAccelRun(configuration, param, self.rundir,
                               parallel=self.do_parallel)
        run.run()

        # Use loaded run
        self.run = maccel.MAccelRun.load(pjoin(run.rundir, 'run.pickl'))
        movie = maccel.PlotMaker(run, load_dir=self.rundir,
                                 parallel=self.do_parallel)
        movie.make_plots()

    def test_files(self):
        nrounds = self.run.n_rounds
        for round_i in range(nrounds * 2):
            should_exist = round_i < nrounds
            with self.subTest(round_i=round_i, should_exist=should_exist):
                plotfn = 'plot-{round_i:04d}.png'.format(round_i=round_i)
                plotfn = pjoin(self.rundir, 'figs', plotfn)
                if should_exist:
                    self.assertTrue(os.path.exists(plotfn))
                else:
                    self.assertFalse(os.path.exists(plotfn))


    def tearDown(self):
        shutil.rmtree(self.rundir)


class TestPlotNoParallel(TestPlot):
    def setUp(self):
        self.do_parallel = False
        super().setUp()


class TestConvergence(TestCase):
    def setUp(self):
        configuration = maccel.AlanineConfiguration().apply_configuration()
        self.tpr = 3
        self.spt = 4000
        param = maccel.AlanineParams(spt=self.spt, tpr=self.tpr, step_res=-1,
                                     post_converge=2)
        self.rundir = get_folder('plot')
        run = maccel.MAccelRun(configuration, param, self.rundir,
                               parallel=False)
        run.run()

        self.runfn = pjoin(run.rundir, 'run.pickl')

    def test_find_convergence_from_fn(self):
        dd = maccel.find_convergence_from_filename(self.runfn)
        self.assertEqual(dd['run_id'], 0)
        self.assertEqual(dd['spt'], self.spt)
        self.assertEqual(dd['tpr'], self.tpr)
        self.assertTrue(dd['steps'] > 0)
        self.assertTrue(dd['rounds'] > 0)


if __name__ == "__main__":
    unittest.main()