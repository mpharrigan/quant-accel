__author__ = 'harrigan'

import os
from os.path import join as pjoin
import logging as log

from IPython.parallel import Client

from .run import MAccelRun
import itertools


class MAccelGrid(object):
    """Run many parameter configurations."""

    def __init__(self, configuration, griddir):
        self.config = configuration
        self.griddir = griddir

        try:
            c = Client()
            self.lbv = c.load_balanced_view()
            self.lbv.block = True
        except FileNotFoundError as e:
            log.error("Could not connect to parallel engine: %s", e)
            self.lbv = None

    def grid(self):
        """Launch several runs over a grid of parameters."""

        if self.lbv is None:
            return False

        self.lbv.map(_launch, zip(itertools.repeat(self.config),
                                  itertools.repeat(self.griddir),
                                  self.config.get_param_grid()))


def _launch(arg_tuple):
    """Helper function for map."""
    config, griddir, params = arg_tuple
    rundir = pjoin(griddir, params.dirname)
    run = MAccelRun(config, params, rundir)
    run.run()
