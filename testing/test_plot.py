__author__ = 'harrigan'
from unittest import TestCase
from os.path import join as pjoin
import os
import logging

import maccelerator as maccel
from maccelerator.test_utils import get_folder, get_fn

# Disable logging during test
logging.captureWarnings(True)
logging.disable(logging.WARNING)


class TestPlot(TestCase):
    def setUp(self):
        configuration = maccel.AlanineConfiguration(get_fn('ala.msm.pickl'),
                                                    get_fn('ala.centers.h5'))
        self.tpr = 3
        self.spt = 1000
        param = maccel.AlanineParams(spt=self.spt, tpr=self.tpr)
        self.rundir = get_folder('plot')
        run = maccel.MAccelRun(configuration, param, self.rundir)
        run.run()

        # Use loaded run
        self.run = maccel.MAccelRun.load(pjoin(run.rundir, 'run.pickl'))
        movie = maccel.PlotMaker(run)
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


