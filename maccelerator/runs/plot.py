"""Code for making plots from a completed run."""
__author__ = 'harrigan'

import logging
from os.path import join as pjoin
import os
import itertools
import pickle

import pandas as pd
import numpy as np
from IPython.parallel import Client

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
        args = [self._get_for_parallel(i) for i in range(n_rounds)]
        self.lbv.map(_plot_helper, args)

    def _get_for_parallel(self, round_i, rel=True):
        """Create a tuple of arguments for parallel helper."""
        file = self.run.config.file
        out_fn = file.plot_fn(round_i, rel=rel)
        out_fn = pjoin(self.load_dir, out_fn)
        out_fn = os.path.abspath(out_fn)

        return self.load_convergence(round_i)[-1], self.run.params, out_fn

    def load_convergence(self, round_i, rel=True):
        """Load a convergence object for a particular round."""
        file = self.run.config.file
        conv_fn = "{}.pickl".format(file.conv_fn(round_i, rel=rel))
        conv_fn = pjoin(self.load_dir, conv_fn)

        with open(conv_fn, 'rb') as conv_f:
            converge = pickle.load(conv_f)
        return converge

    def load_convergences(self):
        """Load all convergences"""
        return [self.load_convergence(i) for i in range(self.run.n_rounds)]

    def convergence_dataframe(self):
        """Get a dataframe of convergences over time."""
        round_is = range(self.run.n_rounds)
        substeps = self.run.params.subbuild_uptos
        coords = np.array(list(itertools.product(round_is, substeps)))
        steps = self.run.params.spt * coords[:, 0] + coords[:, 1]

        conv_vals = np.asarray(
            [[c.converged for c in cs] for cs in self.load_convergences()]
        ).reshape(-1)

        df = pd.DataFrame(dict(
            round_i=coords[:, 0], steps=steps, converged=conv_vals
        )).set_index('steps')

        return df

    def find_first_convergence(self, window=4, cutoff=0.5):
        """Use a rolling average to find step and round of first convergence.
        """
        conv_df = self.convergence_dataframe()
        rolling_df = pd.rolling_mean(conv_df['converged'], window).fillna(0)
        steps = (rolling_df >= cutoff).argmax()
        rounds = conv_df['round_i'].loc[steps] + 1
        return steps, rounds


def _plot_helper(args):
    """Can be mapped."""
    converge, params, fn = args
    converge.plot_and_save(params, None, fn)


def find_convergence_from_filename(run_fn):
    """From a filename, return convergence data suitable for pandas

    Use this from IPython.parallel map
    """

    with open(run_fn, 'rb') as run_f:
        run = pickle.load(run_f)
    pm = PlotMaker(run, load_dir=os.path.dirname(run_fn), parallel=False)
    steps, rounds = pm.find_first_convergence()

    return dict(run_fn=run_fn, run_id=run.params.run_id,
                spt=run.params.spt, tpr=run.params.tpr,
                steps=steps, rounds=rounds)