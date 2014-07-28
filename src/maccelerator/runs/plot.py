"""Code for making plots from a completed run."""
__author__ = 'harrigan'

import logging

from IPython.parallel import Client

from ..convergence.base import Convergence
from .run import NoParallelView


log = logging.getLogger(__name__)


class PlotMaker():
    """From a completed run, produce plots.

    :param run: The run object
    :param parallel: Whether to use IPython parallel api
    """

    def __init__(self, run, parallel=True):
        self.run = run
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

        log.debug('Making %d frames', n_rounds)
        self.lbv.map(_plot_helper,
                     [self._get_for_parallel(i) for i in range(n_rounds)])

    def _get_for_parallel(self, round_i):
        """Create a tuple of arguments for parallel helper."""
        file = self.run.config.file
        converge = Convergence.load("{}.pickl".format(file.conv_fn(round_i)))
        return converge, self.run.params, file.plot_fn(round_i)


def _plot_helper(args):
    """Can be mapped."""
    converge, params, fn = args
    converge.plot_and_save(params, None, fn)
