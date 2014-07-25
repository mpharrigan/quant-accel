"""
Classes for building models

Created on Mar 5, 2014

@author: harrigan
"""

import logging
import pickle

import numpy as np
from mixtape.cluster import MiniBatchKMeans
from mixtape.markovstatemodel import MarkovStateModel

log = logging.getLogger(__name__)


class Modeller:
    """Base class for constructing models."""

    def model(self, traj_fns, params):
        """Create a model.

        :param traj_fns: Trajectory filenames from which we build the model
        :param params: Simulation parameters.
        """
        raise NotImplementedError

    @property
    def modelfn(self):
        return "msm-{round_i}"


class Model:
    def __init__(self, msm):
        self._msm = msm

    def save(self, fn):
        fn = "{}.pickl".format(fn)
        with open(fn, 'wb') as f:
            pickle.dump(self, f)

    @property
    def n_states(self):
        return self._msm.n_states_

    @property
    def counts(self):
        return self._msm.raw_counts_

    @property
    def timescales(self):
        return self._msm.timescales_

    @property
    def lagtime(self):
        return self._msm.lag_time


class ClusterModel(Model):
    def __init__(self, msm, clusterer):
        super().__init__(msm)
        self._clusterer = clusterer


class ClusterModeller(Modeller):
    """Cluster and then build a model from cluster labels."""


    def load_trajs(self, traj_fns):
        raise NotImplementedError

    def lagtime(self, params):
        return params.build_lt


    def model(self, traj_fns, params):
        """Cluster using kmeans and build an MSM


        References:
        http://en.wikipedia.org/wiki
        /Determining_the_number_of_clusters_in_a_data_set#Rule_of_thumb
        """
        log.debug("Loading trajectories")
        trajs_vec = self.load_trajs(traj_fns)

        log.info("Starting cluster")

        # Get number of data points
        n_datapoints = np.sum([len(traj) for traj in trajs_vec])
        n_clusters = int(np.sqrt(n_datapoints / 2))

        # Do clustering
        clusterer = MiniBatchKMeans(n_clusters=n_clusters)
        clusterer.fit(trajs_vec)

        log.info("Building MSM")
        msm = MarkovStateModel(lag_time=self.lagtime(params), n_timescales=10)
        msm.fit(clusterer.labels_)

        return ClusterModel(msm, clusterer)


class TMatModel(Model):
    def __init__(self, msm, full_tmat, full_populations, full_eigenvec,
                 found_states):
        super().__init__(msm)
        self.full_tmat = full_tmat
        self.full_populations = full_populations
        self.full_eigenvec = full_eigenvec
        self.found_states = found_states
        self.tot_n_states = len(full_populations)

    @property
    def timescales(self):
        if self._msm.n_states_ >= 2:
            return self._msm.timescales_
        else:
            return np.asarray([0])


class TMatModeller(Modeller):
    """Model from transition matrix trajectories. (No clustering)"""

    def load_trajs(self, traj_fns):
        raise NotImplementedError

    def lagtime(self, params):
        raise NotImplementedError

    def __init__(self, tot_n_states):
        super().__init__()
        self.tot_n_states = tot_n_states

    def model(self, traj_fns, params):
        """Build a model from the result of a transition matrix simulations

        We take care of only returning states which we have 'discovered'

        """

        msm = MarkovStateModel(lag_time=self.lagtime(params))
        msm.fit(self.load_trajs(traj_fns))

        log.debug("Number of untrimmed States: %d", msm.n_states_)

        # Back out full-sized populations
        if msm.n_states_ < self.tot_n_states:
            populations = np.zeros(self.tot_n_states)
            eigenvec = np.zeros(self.tot_n_states)
            tmat = np.zeros((self.tot_n_states, self.tot_n_states))

            kept = np.array(msm.state_labels_)

            # Be super explicit about minimum state requirements
            if msm.n_states_ >= 1:
                populations[kept] = msm.populations_
                tmat[np.ix_(kept, kept)] = msm.transmat_

            if msm.n_states_ >= 2:
                eigenvec[kept] = msm.left_eigenvectors_[:, 1]

        else:
            populations = msm.populations_
            eigenvec = msm.left_eigenvectors_[:, 1]
            tmat = msm.transmat_

        # Get found states
        # Those which have at least one transition to or from
        # Note: We can't just sample from states with zero 'from' counts
        # This would neglect states visited at the ends of trajectories.
        # These are probably pretty important for adaptive sampling
        # countscoo = msm.rawcounts_.tocoo()
        # found_states = np.hstack((countscoo.row, countscoo.col))
        # found_states = np.unique(found_states)
        # TODO
        found_states = None

        return TMatModel(msm, tmat, populations, eigenvec, found_states)

