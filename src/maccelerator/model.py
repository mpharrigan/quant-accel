"""
Classes for building models and adapting from them

Created on Mar 5, 2014

@author: harrigan
"""

import logging as log

import numpy as np
from mixtape.cluster import MiniBatchKMeans
from mixtape.markovstatemodel import MarkovStateModel


class Adapter(object):
    """Base class for an object that chooses where to start new simulation."""

    def adapt(self, params):
        """Return a state from which to start.

        :param params: Simulation parameters.
        """
        raise NotImplementedError


class SortCountsAdapter(Adapter):
    """Choose the states from which we've transitioned the fewest times,
    in order
    """

    def __init__(self, modeller):
        super().__init__()
        self.modeller = modeller

    def adapt(self, params):
        """From a counts matrix, pick the best states from which to start.

        :param params: Simulation parameters so we know how many new states
                       to return
        :returns: Index of state
        """
        counts = self.modeller.counts
        found_states = None  # TODO: Deal with found states
        n_tpr = params.tpr

        # Get counts
        counts_per_state = np.asarray(counts.sum(axis=1)).flatten()

        # Only consider states we know about
        if found_states is not None:
            counts_per_state = counts_per_state[found_states]

        # Sort
        states_to_sample = np.argsort(counts_per_state)

        # Get the right number of states
        if len(states_to_sample) > n_tpr:
            states_to_sample = states_to_sample[:n_tpr]

        log.info('Generating %d new starting structures.',
                 len(states_to_sample))
        counts_str = ', '.join(
            [str(j) for j in counts_per_state[states_to_sample]])
        log.debug('Counts %s', counts_str)
        return states_to_sample


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

    def __init__(self):
        super().__init__()
        self.msm = None
        self.found_states = None

    def _model(self, trajs, lagtime):
        """Build a model from the result of a transition matrix simulations

        We take care of only returning states which we have 'discovered'

        :param trajs: List of ndarray; State indices
        :param lagtime: Build a model at this lag time
        """

        msm = MarkovStateModel(lag_time=lagtime, n_timescales=10)
        msm.fit(trajs)
        self.msm = msm

        # Get found states
        # Those which have at least one transition to or from
        # Note: We can't just sample from states with zero 'from' counts
        # This would neglect states visited at the ends of trajectories.
        # These are probably pretty important for adaptive sampling
        countscoo = msm.rawcounts_.tocoo()
        found_states = np.hstack((countscoo.row, countscoo.col))
        self.found_states = np.unique(found_states)

    @property
    def counts(self):
        """Raw counts to use for adapting."""
        return self.msm.rawcounts_
