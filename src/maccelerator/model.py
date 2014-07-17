"""
Created on Mar 5, 2014

@author: harrigan
"""

import logging as log
import os

import mdtraj as md
import numpy as np

from mdtraj import io

from mixtape.cluster import KMeans, MiniBatchKMeans
from mixtape.markovstatemodel import MarkovStateModel


class Adapter(object):
    def adapt(self, params):
        raise NotImplementedError


class SortCountsAdapter(Adapter):
    def __init__(self, modeller):
        super().__init__()
        self.modeller = modeller

    def adapt(self, params):
        """From a counts matrix, pick the best states from which to start.

        :returns: Index of state
        """
        counts = self.modeller.counts
        found_states = None  #TODO
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
    def model(self, traj_fns, params):
        raise NotImplementedError

    def seed_state(self, params):
        """Get seed states to start the run.

        :param params: Contains number of seed states to generate
        """
        raise NotImplementedError


class ClusterModeller(Modeller):
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
        clusterer = MiniBatchKMeans(n_clusters=n_clusters)
        clusterer.fit(trajs_vec)
        self.clusterer = clusterer

        log.info("Building MSM")
        msm = MarkovStateModel(lag_time=lagtime, n_timescales=10)
        msm.fit(clusterer.labels_)
        self.msm = msm

    @property
    def counts(self):
        return self.msm.rawcounts_


class TMatModeller(Modeller):
    def tmat_model(self, trajs, lagtime):
        """Get counts from a tmat simulation

        We take care of only returning states which we have 'discovered'

        :param trajs: List of ndarray
        :param lagtime: Lagtime
        :return: counts: counts of states we have discovered.
                 found_states: actual indices of the states, used for translating
                 back into absolute state indices
        """

        # Don't need n_states, we won't be sampling from anything higher than
        # that anyways
        counts = msml.get_count_matrix_from_assignments(np.array(trajs),
                                                        lag_time=lagtime,
                                                        sliding_window=True)

        # Get found states
        # Those which have at least one transition to or from
        # Note: We can't just sample from states with zero 'from' counts
        # This would neglect states visited at the ends of trajectories.
        # These are probably pretty important for adaptive sampling
        countscoo = counts.tocoo()
        found_states = np.hstack((countscoo.row, countscoo.col))
        found_states = np.unique(found_states)
        return counts, found_states


# TODO: Move these
class MoveTheseElsewhere(object):
    def load_muller(in_fn):
        """Use mdtraj for loading these trajectories."""
        return md.load(in_fn)


    def save_muller(out_fn, traj):
        """Save a trajectory of centroids."""
        traj.save(out_fn)


    def load_tmat(in_fn):
        """Use io.loadh to load tmat state indices."""
        return io.loadh(in_fn, 'arr_0')


    def save_tmat(out_fn, state_is):
        """Save a matrix."""
        io.saveh(out_fn, state_is)


    def save_starting_states(state_is, round_i, save_func):
        """Save to a consistent filename.

        :param save_func: Use this function to actually save.
        """
        out_fn = os.path.join('sstates', 'round-%d.h5' % (round_i + 1))
        save_func(out_fn, state_is)


    def load_trajectories(round_i, load_func):
        """Load trajectories up to and including :round_i:

        Helper function for model building.
        """

        trajs = []
        for cround in range(round_i + 1):
            tdir = os.path.join('trajs', 'round-%d' % cround)
            trajs += [load_func(os.path.join(tdir, s))
                      for s in os.listdir(tdir) if s.endswith('.h5')]

        # Stats
        # Note: Use len(), which works on mdtraj trajectories and ndarrays
        traj_len = len(trajs[0])
        wall_steps = traj_len * (round_i + 1)
        return wall_steps, trajs


    def load_trajectories_percent(percent, load_func):
        """Load trajectories up to a certain percent.

        Helper function for model building, esp. for dealing with long
        trajectory (non-adaptive) runs against which we wish to compare.
        """

        trajs = []
        tdir = os.path.join('trajs', 'round-0')
        trajs += [load_func(os.path.join(tdir, s))
                  for s in os.listdir(tdir) if s.endswith('.h5')]

        # Find ending index
        # Note: Use len(), which works on md.Trajectories and ndarrays
        traj_len = len(trajs[0])
        endstep = int(percent * traj_len)

        trajs = [traj[:endstep] for traj in trajs]
        return endstep, trajs


    def model_and_adapt_tmat(args):
        """Load trajectories, model, and select starting states for tmat.

        :param args: Arguments from argparse.
        """
        _, trajs = load_trajectories(args.round, load_tmat)
        counts, found_states = tmat_model(trajs, args.lagtime)
        sstates_sub = adapt(counts, args.n_tpr, found_states)
        # Translate indices from found_states into absolute indices
        sstates = found_states[sstates_sub]
        save_starting_states(sstates, args.round, save_tmat)


    def model_and_adapt_muller(args):
        """Load trajectories, model, and select starting states for muller.

        :param args: Arguments from argparse
        """
        _, trajs = load_trajectories(args.round, load_muller)
        counts, centroids = cluster_model(trajs, args.lagtime,
                                          args.distance_cutoff)
        sstates = adapt(counts, args.n_tpr)
        # Translate state indices to centroids
        traj = centroids[sstates]
        save_starting_states(traj, args.round, save_muller)
