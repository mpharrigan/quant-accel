__author__ = 'harrigan'

from os.path import join as pjoin
import itertools
import os
import logging
import glob
import shutil
import traceback

from IPython.parallel import Client

from .run import MAccelRun


log = logging.getLogger(__name__)
SHM = '/dev/shm/'


class MAccelGrid:
    """Run many parameter configurations.

    :param parallel: [True, False]
        Whether to run each *run* as a parallel thing.
        Support for running grids in parallel has been dropped.
    """

    def __init__(self, configuration, griddir, run_id=0, parallel=False):
        self.config = configuration
        self.griddir = griddir
        self.run_id = run_id
        self.grid = self._grid_noparallel

        self.run_parallel = parallel

    def _grid_parallel(self):
        """Launch several runs over a grid of parameters."""
        log.warning("_grid_parallel is deprecated!")

        # Check
        if self.lbv is None:
            return False

        # Run helper function
        self.lbv.map(_launch, zip(itertools.repeat(self.config),
                                  itertools.repeat(self.griddir),
                                  self.config.get_param_grid(self.run_id),
                                  itertools.repeat(self.run_parallel)))

    def _grid_noparallel(self):
        """Do runs serially"""
        for arg_tuple in zip(itertools.repeat(self.config),
                             itertools.repeat(self.griddir),
                             self.config.get_param_grid(self.run_id),
                             itertools.repeat(self.run_parallel)):
            try:
                _launch(arg_tuple)
            except Exception as e:
                log.error(str(e))
                log.error(traceback.format_exc())


def _launch(arg_tuple):
    """Helper function for map."""
    config, griddir, params, parallel = arg_tuple
    rundir = pjoin(griddir, params.dirname)

    # Make the directory for the run
    try:
        os.mkdir(rundir)
    except OSError as e:
        log.error('Problem with %s: %s; skipping!', rundir, e)

    log.info('Launching run %s', params.pretty_desc)
    run = MAccelRun(config, params, rundir, parallel)
    run.run()


def _archive_trajs(griddir):
    for traj_dir in glob.iglob(pjoin(griddir, '*/trajs/')):
        logging.info('Archiving %s', traj_dir)
        shutil.make_archive(base_name=pjoin(traj_dir, '..', 'trajs'),
                            format='gztar',
                            root_dir=pjoin(traj_dir, '..'),
                            base_dir='trajs')
        shutil.rmtree(traj_dir)


class MaccelGridFs:
    """Use filesystem."""

    def __init__(self, config, run_id):
        self.config = config
        self.run_id = run_id
        self.grid = None

    def __enter__(self):
        self.griddir = 'copy-{}'.format(self.run_id)

        assert not os.path.exists(self.griddir)
        os.mkdir(self.griddir)

        self.grid = MAccelGrid(self.config, self.griddir, self.run_id,
                               parallel=False)
        return self.grid

    def __exit__(self, exc_type, exc_val, exc_tb):
        _archive_trajs(self.griddir)


class MaccelGridShm:
    """Use /dev/shm/ to reduce disk io."""

    def __init__(self, config, run_id):
        self.config = config
        self.run_id = run_id

        self.grid = None
        self.griddir = ""
        self.shm_griddir = ""

    def __enter__(self):
        self.griddir = 'copy-{}'.format(self.run_id)
        self.shm_griddir = pjoin(SHM, self.griddir)

        # Do it in shared memory
        assert not os.path.exists(self.griddir)
        os.mkdir(self.shm_griddir)

        # Initialize
        self.grid = MAccelGrid(self.config, self.shm_griddir, self.run_id,
                               parallel=False)
        return self.grid


    def __exit__(self, exc_type, exc_val, exc_tb):
        _archive_trajs(self.shm_griddir)
        shutil.copytree(self.shm_griddir, self.griddir)
        shutil.rmtree(self.shm_griddir)
