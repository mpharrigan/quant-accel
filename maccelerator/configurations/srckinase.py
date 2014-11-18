"""Src Kinase Transition matrix system

Author: Matthew Harrigan
"""

import logging
import itertools

import scipy.io
from mixtape.markovstatemodel import MarkovStateModel
import mdtraj as md
from mixtape.pca import PCA
from mixtape.featurizer import DihedralFeaturizer
import numpy as np

from ..simulate import TMatSimulator
from ..model import TMatModeller
from ..adapt import RandomAdapter
from .base import TMatConfiguration
from ..convergence import TMatConvergenceChecker
from ..param import AdaptiveParams
from ..files import get_fn


log = logging.getLogger(__name__)


def generate_srckinase_msm(tmat_fn, populations_fn, mapping_fn, gens_fn):
    msm = _generate_msm(tmat_fn, populations_fn)
    centers = _generate_centers(mapping_fn, gens_fn)
    return msm, centers


def _generate_msm(tmat_fn, populations_fn):
    tmat_sparse = scipy.io.mmread(tmat_fn)
    tmat_dense = tmat_sparse.toarray()

    populations = np.loadtxt(populations_fn)

    # Data indicates a lag time of 5. Note this is not really relevant
    # Make sure in convergence checks to divide timescales by this
    # for an accurate comparison
    msm = MarkovStateModel(lag_time=5)
    msm.n_states_ = tmat_dense.shape[0]
    msm.mapping_ = dict(zip(np.arange(msm.n_states_), np.arange(msm.n_states_)))
    msm.transmat_ = tmat_dense
    msm.populations_ = populations

    # Force eigensolve and check consistency
    computed_pops = msm.left_eigenvectors_[:, 0]
    computed_pops /= np.sum(computed_pops)

    np.testing.assert_allclose(computed_pops, msm.populations_)

    return msm


def _generate_centers(mapping_fn, gens_fn):
    mapping = np.loadtxt(mapping_fn)

    gens = md.load(gens_fn)
    gens = gens[mapping != -1]

    dihed = DihedralFeaturizer(['phi', 'psi'])
    dihedx = dihed.fit_transform([gens])

    pca = PCA(n_components=2)
    pcax = pca.fit_transform(dihedx)[0]

    return pcax


class SrcKinaseParams(AdaptiveParams):
    """Returned by get_param_grid which must be implemented."""

    def __init__(self, spt, tpr, adapt_lt=1, build_lt=1,
                 post_converge=10, run_id=0, step_res=8):
        super().__init__(spt=spt, tpr=tpr, adapt_lt=adapt_lt, build_lt=build_lt,
                         post_converge=post_converge, run_id=run_id,
                         step_res=step_res)


class SrcKinaseConfiguration(TMatConfiguration):
    def __init__(self):
        super().__init__(get_fn('src.msm.pickl'), get_fn('src.centers.h5'))

    def defaults(config):
        # Simulator
        config.simulator_class = TMatSimulator

        # Modeller
        config.modeller_class = TMatModeller

        # Convergence checker
        config.convchecker_class = TMatConvergenceChecker

        # Adapter
        config.adapter_class = RandomAdapter

        # Define a function to yield combinations of parameters
        def get_param_grid(run_id):
            spts = [10, 100]
            tprs = [1, 10]

            for spt, tpr in itertools.product(spts, tprs):
                yield SrcKinaseParams(spt=spt, tpr=tpr, run_id=run_id)

        config.get_param_grid = get_param_grid

        # Tolerance scaling factor
        config.tolerance_scale = 1.0


