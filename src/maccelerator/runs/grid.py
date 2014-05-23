__author__ = 'harrigan'

import os
from os.path import join as pjoin
import logging as log

from .run import MAccelRun


class MAccelGrid(object):
    """Run many parameter configurations."""

    def __init__(self, configuration, griddir):
        self.config = configuration
        self.griddir = griddir

    def grid(self):
        """Launch several runs over a grid of parameters."""

        for params in self.config.get_param_grid():
            rundir = pjoin(self.griddir, params.dirname)
            try:
                os.mkdir(rundir)
            except OSError as e:
                log.warning("Skipping %s (%s)", rundir, e)
                continue

            run = MAccelRun(self.config, params, rundir)
            run.run()
