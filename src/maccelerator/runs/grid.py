__author__ = 'harrigan'

from os.path import join as pjoin
import logging as log
import itertools
import os

from IPython.parallel import Client

from .run import MAccelRun


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

        # Check
        if self.lbv is None:
            return False


        # Run helper function
        self.lbv.map(_launch, zip(itertools.repeat(self.config),
                                  itertools.repeat(self.griddir),
                                  self.config.get_param_grid()))

    def grid_noparallel(self):
        """Do runs serially (for debugging)."""
        for arg_tuple in zip(itertools.repeat(self.config),
                             itertools.repeat(self.griddir),
                             self.config.get_param_grid()):
            _launch(arg_tuple)


def _launch(arg_tuple):
    """Helper function for map."""
    config, griddir, params = arg_tuple
    rundir = pjoin(griddir, params.dirname)

    # Make the directory for the run
    try:
        os.mkdir(rundir)
    except OSError:
        log.error('%s Already exists, skipping!', rundir)

    run = MAccelRun(config, params, rundir)
    run.run()
