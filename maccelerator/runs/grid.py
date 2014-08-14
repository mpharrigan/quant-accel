__author__ = 'harrigan'

from os.path import join as pjoin
import itertools
import os
import logging

from IPython.parallel import Client

from .run import MAccelRun


log = logging.getLogger(__name__)


class MAccelGrid(object):
    """Run many parameter configurations.

    :param parallel: [True, False] or ['full', 'semi', 'non']
        'semi' will do parallel simulation but not parallel over grid
    """

    def __init__(self, configuration, griddir, run_id=0, parallel=True):
        self.config = configuration
        self.griddir = griddir
        self.run_id = run_id

        if parallel or parallel == 'full':
            self.grid = self._grid_parallel
            try:
                c = Client()
                self.lbv = c.load_balanced_view()
                self.lbv.block = True
            except FileNotFoundError as e:
                log.error("Could not connect to parallel engine: %s", e)
                self.lbv = None
        else:
            self.grid = self._grid_noparallel

        if not parallel or parallel == 'non':
            self.run_parallel = False
        else:
            self.run_parallel = True

    def _grid_parallel(self):
        """Launch several runs over a grid of parameters."""

        # Check
        if self.lbv is None:
            return False

        # Run helper function
        self.lbv.map(_launch, zip(itertools.repeat(self.config),
                                  itertools.repeat(self.griddir),
                                  self.config.get_param_grid(self.run_id),
                                  itertools.repeat(self.run_parallel)))

    def _grid_noparallel(self):
        """Do runs serially (for debugging)."""
        for arg_tuple in zip(itertools.repeat(self.config),
                             itertools.repeat(self.griddir),
                             self.config.get_param_grid(self.run_id),
                             itertools.repeat(self.run_parallel)):
            _launch(arg_tuple)


def _launch(arg_tuple):
    """Helper function for map."""
    config, griddir, params, parallel = arg_tuple
    rundir = pjoin(griddir, params.dirname)

    # Make the directory for the run
    try:
        os.mkdir(rundir)
    except OSError:
        log.error('%s Already exists, skipping!', rundir)

    log.info('Launching run %s', params.pretty_desc)
    run = MAccelRun(config, params, rundir, parallel)
    run.run()
