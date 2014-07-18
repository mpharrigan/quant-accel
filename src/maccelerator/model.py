"""
Classes for building models

Created on Mar 5, 2014

@author: harrigan
"""

import logging as log

import numpy as np
from mixtape.cluster import MiniBatchKMeans
from mixtape.markovstatemodel import MarkovStateModel
import scipy.sparse

# Minimum number of states to have succeeded in building a model
MINSTATES = 4


class Modeller(object):
    """Base class for constructing models."""

    def model(self, traj_fns, params):
        """Create a model.

        :param traj_fns: Trajectory filenames from which we build the model
        :param params: Simulation parameters.
        """
        raise NotImplementedError

    def seed_state(self, params):
        """Get seed states to start the run.

        :param params: Contains number of seed states to generate
        """
        raise NotImplementedError


class ClusterModeller(Modeller):
    """Cluster and then build a model from cluster labels."""

    def __init__(self):
        super().__init__()
        self.msm = None
        self.clusterer = None

    def _model(self, trajs_vec, lagtime):
        """Cluster using kmeans and build an MSM

        :param trajs_vec: List of vector representation of trajectories.
        :param lagtime: The desired lagtime of the model.

        References:
        http://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set#Rule_of_thumb
        """
        log.info("Starting cluster")

        # Get number of data points
        n_datapoints = np.sum([len(traj) for traj in trajs_vec])
        n_clusters = int(np.sqrt(n_datapoints / 2))

        # Do clustering
        clusterer = MiniBatchKMeans(n_clusters=n_clusters)
        clusterer.fit(trajs_vec)
        self.clusterer = clusterer

        log.info("Building MSM")
        msm = MarkovStateModel(lag_time=lagtime, n_timescales=10)
        msm.fit(clusterer.labels_)
        self.msm = msm

    @property
    def counts(self):
        """Raw counts from which we can estimate uncertainty for adapting."""
        return self.msm.rawcounts_


class TMatModeller(Modeller):
    """Model from transition matrix trajectories. (No clustering)"""

    def __init__(self, tot_n_states):
        super().__init__()
        self.msm = None
        self.found_states = None
        self.full_populations = None
        self.full_eigenvec = None
        self.full_tmat = None
        self.tot_n_states = tot_n_states

    def _model(self, trajs, lagtime):
        """Build a model from the result of a transition matrix simulations

        We take care of only returning states which we have 'discovered'

        :param trajs: List of ndarray; State indices
        :param lagtime: Build a model at this lag time
        """

        msm = MarkovStateModel(lag_time=lagtime,
                               n_states=self.tot_n_states)
        msm.fit(trajs)
        self.msm = msm

        log.debug("States left: %d", msm.transmat_.shape[0])
        if msm.transmat_.shape[0] < MINSTATES:
            return False

        # Back out full-sized populations
        n_states = self.msm.n_states
        if msm.transmat_.shape[0] < n_states:
            populations = np.zeros(n_states)
            eigenvec = np.zeros(n_states)
            for i in range(n_states):
                try:
                    populations[i] = msm.populations_[msm.mapping_[i]]
                    eigenvec[i] = msm.eigenvectors_[msm.mapping_[i], 1]
                except KeyError:
                    pass
        else:
            populations = msm.populations_
            eigenvec = msm.eigenvectors_[:, 1]

        # Back out full transition matrix, oh boy
        tmatcoo = msm.transmat_.tocoo()
        for fr, to in msm.mapping_.items():
            tmatcoo.row[tmatcoo.row == to] = fr
            tmatcoo.col[tmatcoo.col == to] = fr

        self.full_populations = populations
        self.full_eigenvec = eigenvec
        self.full_tmat = scipy.sparse.coo_matrix(
            (tmatcoo.data, (tmatcoo.row, tmatcoo.col))).tocsr()

        # Get found states
        # Those which have at least one transition to or from
        # Note: We can't just sample from states with zero 'from' counts
        # This would neglect states visited at the ends of trajectories.
        # These are probably pretty important for adaptive sampling
        countscoo = msm.rawcounts_.tocoo()
        found_states = np.hstack((countscoo.row, countscoo.col))
        self.found_states = np.unique(found_states)

        return True

    @property
    def counts(self):
        """Raw counts to use for adapting."""
        return self.msm.rawcounts_
