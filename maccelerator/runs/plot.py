"""Code for making plots from a completed run."""
__author__ = 'harrigan'

import logging
from os.path import join as pjoin

from IPython.parallel import Client

from ..convergence.base import SupConvergence
from .run import NoParallelView


log = logging.getLogger(__name__)


class PlotMaker():
    """From a completed run, produce plots.

    :param run: The run object
    :param parallel: Whether to use IPython parallel api
    """

    def __init__(self, run, parallel=True, load_dir='.'):
        self.run = run
        self.load_dir = load_dir

        # Set up optional parallelization
        if parallel:
            try:
                c = Client()
                self.lbv = c.load_balanced_view()
                self.lbv.block = True
            except FileNotFoundError as e:
                log.error("Could not connect to parallel engine: %s", e)
                self.lbv = None
        else:
            self.lbv = NoParallelView()

    def make_plots(self):
        """Make plots for all rounds."""
        n_rounds = self.run.n_rounds

        log.info('Making %d frames', n_rounds)
        self.lbv.map(_plot_helper,
                     [self._get_for_parallel(i, rel=True) for i in
                      range(n_rounds)])

    def _get_for_parallel(self, round_i, rel=False):
        """Create a tuple of arguments for parallel helper."""

        file = self.run.config.file
        out_fn = file.plot_fn(round_i, rel=rel)
        conv_fn = "{}.pickl".format(file.conv_fn(round_i, rel=rel))

        if rel:
            conv_fn = pjoin(self.load_dir, conv_fn)
            out_fn = pjoin(self.load_dir, out_fn)

        converge = SupConvergence.load(conv_fn)
        return converge, self.run.params, out_fn

    def load_convergences(self):
        return [self._get_for_parallel(i, rel=True)[0] for i in
                range(self.run.n_rounds)]


def _plot_helper(args):
    """Can be mapped."""
    converge, params, fn = args
    converge.plot_and_save(params, None, fn)
