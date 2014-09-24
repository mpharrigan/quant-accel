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
import mdtraj
import scipy.sparse


log = logging.getLogger(__name__)


class Modeller:
    """Base class for constructing models."""

    def __init__(self, config):
        # pycharm: disable=unused-argument
        self.n_builds = 1

    def multi_model(self, traj_fns, params, step_res=None):
        """Build several models for one round."""
        n_builds = params.spt // step_res
        assert params.spt % step_res == 0, 'Must divide evenly'
        up_tos = (np.arange(n_builds) + 1) * step_res
        return [self.model(traj_fns, params, up_to=up_to) for up_to in up_tos]

    def model(self, traj_fns, params, up_to=None):
        """Create a model.

        :param traj_fns: Trajectory filenames from which we build the model
        :param params: Simulation parameters.
        """
        raise NotImplementedError

    @property
    def modelfn(self):
        return "msm-{round_i}"


class Model:
    def __init__(self, params):
        self._n_states = 0
        self._tmat = np.zeros((0, 0))
        self._populations = np.zeros(0)
        self._timescales = np.zeros(0)
        self._eigenvectors = np.zeros((0, 0))
        self._eigenvalues = np.zeros(0)
        self._adapt_counts = np.zeros((0, 0))
        self._full_counts = np.zeros((0, 0))

        self.params = params

        self._is_consistent = True
        self._dirty = False


    def save(self, fn):
        fn = "{}.pickl".format(fn)
        with open(fn, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, fn):
        with open(fn, 'rb') as f:
            return pickle.load(f)

    def _debug(self):
        newd = dict()
        for k, v in self.__dict__.items():
            if k.startswith('_'):
                try:
                    newd[k] = len(v)
                except:
                    pass
        return newd


    def _check_consistency(self):
        """Check that each item has the same number of states

        Ignored: number of columns in eigenvectors
                 adapt_counts
        """
        if self._is_consistent and not self._dirty:
            return
        elif not self._is_consistent and not self._dirty:
            raise AssertionError("Inconsistent number of states")

        debug_msg = "Inconsistent number of states: {}".format(self._debug())
        ns = [self._tmat.shape[0], self._tmat.shape[1], len(self._populations),
              self._eigenvectors.shape[0], self._full_counts.shape[0],
              self._full_counts.shape[1]]
        assert all(n == ns[0] for n in ns), debug_msg
        self._is_consistent = True
        self._dirty = False
        self._n_states = ns[0]

    @property
    def n_states(self):
        self._check_consistency()
        return self._n_states

    @property
    def tmat(self):
        self._check_consistency()
        return self._tmat

    @tmat.setter
    def tmat(self, value):
        self._tmat = value
        self._dirty = True

    @property
    def populations(self):
        self._check_consistency()
        return self._populations

    @populations.setter
    def populations(self, value):
        self._populations = value
        self._dirty = True

    @property
    def timescales(self):
        self._check_consistency()
        return self._timescales

    @timescales.setter
    def timescales(self, value):
        self._timescales = value
        self._dirty = True

    @property
    def eigenvectors(self):
        self._check_consistency()
        return self._eigenvectors

    @eigenvectors.setter
    def eigenvectors(self, value):
        self._eigenvectors = value
        self._dirty = True

    @property
    def eigenvalues(self):
        self._check_consistency()
        return self._eigenvalues

    @eigenvalues.setter
    def eigenvalues(self, value):
        self._eigenvalues = value
        self._dirty = True

    @property
    def adapt_counts(self):
        self._check_consistency()
        return self._adapt_counts

    @adapt_counts.setter
    def adapt_counts(self, value):
        self._adapt_counts = value
        self._dirty = True

    @property
    def full_counts(self):
        self._check_consistency()
        return self._full_counts

    @full_counts.setter
    def full_counts(self, value):
        self._full_counts = value
        self._dirty = True


class ClusterModel(Model):
    # TODO: Implement
    pass


class ClusterModeller(Modeller):
    """Cluster and then build a model from cluster labels."""

    def load_trajs(self, traj_fns):
        raise NotImplementedError

    def lagtime(self, params):
        return params.build_lt


    def model(self, traj_fns, params, up_to=None):
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
        msm = MarkovStateModel(lag_time=self.lagtime(params), n_timescales=10,
                               verbose=log.level > logging.DEBUG)
        msm.fit(clusterer.labels_)

        return ClusterModel(msm, clusterer)


class TMatModel(Model):
    def __init__(self, params):
        super().__init__(params)
        self.found_states = None
        self.adapt_mapping = None


class TMatModeller(Modeller):
    """Model from transition matrix trajectories. (No clustering)"""

    def __init__(self, config):
        super().__init__(config)
        self.tot_n_states = config.ref_msm.n_states_

    def load_trajs(self, traj_fns, up_to=None):
        if up_to is not None:
            trajs = [mdtraj.io.loadh(fn, 'state_traj')[:up_to] for fn in
                     traj_fns]
        else:
            trajs = [mdtraj.io.loadh(fn, 'state_traj') for fn in traj_fns]
        return trajs


    def model(self, traj_fns, params, up_to=None):
        """Build a model from the result of a transition matrix simulations

        We take care of only returning states which we have 'discovered'

        """
        ret_model = TMatModel(params)
        trajs = self.load_trajs(traj_fns, up_to)

        # -----------------------------------------------------------
        # Do building for adaptive
        # -----------------------------------------------------------
        msm_adapt = MarkovStateModel(lag_time=params.adapt_lt, verbose=False,
                                     ergodic_cutoff=0, prior_counts=1,
                                     reversible_type='none')
        msm_adapt.fit(trajs)
        ret_model.adapt_counts = msm_adapt.countsmat_
        ret_model.adapt_mapping = msm_adapt.mapping_

        if msm_adapt.n_states_ < self.tot_n_states:
            counts = np.zeros((self.tot_n_states, self.tot_n_states))
            kept = np.array(msm_adapt.state_labels_)
            counts[np.ix_(kept, kept)] = msm_adapt.countsmat_
            ret_model.full_counts = counts
        else:
            ret_model.full_counts = msm_adapt.countsmat_

        # Get found states
        ret_model.found_states = msm_adapt.state_labels_

        # -----------------------------------------------------------
        # Do building for checking convergence
        # -----------------------------------------------------------
        msm_build = MarkovStateModel(lag_time=params.build_lt, verbose=False)
        msm_build.fit(trajs)

        log.debug("Number of untrimmed States: %d", msm_build.n_states_)

        # Back out full-sized populations
        if msm_build.n_states_ < self.tot_n_states:
            populations = np.zeros(self.tot_n_states)
            eigenvec = np.zeros(self.tot_n_states)
            tmat = np.zeros((self.tot_n_states, self.tot_n_states))

            kept = np.array(msm_build.state_labels_)

            # Be super explicit about minimum state requirements
            if msm_build.n_states_ >= 1:
                populations[kept] = msm_build.populations_
                tmat[np.ix_(kept, kept)] = msm_build.transmat_
            else:
                raise ValueError(
                    "{} states? What is this?".format(msm_build.n_states_))

            if msm_build.n_states_ >= 2:
                eigenvec[kept] = msm_build.left_eigenvectors_[:, 1]

        else:
            populations = msm_build.populations_
            eigenvec = msm_build.left_eigenvectors_[:, 1]
            tmat = msm_build.transmat_

        ret_model.tmat = tmat
        ret_model.populations = populations
        ret_model.timescales = msm_build.timescales_
        ret_model.eigenvectors = eigenvec
        ret_model.eigenvalues = msm_build.eigenvalues_

        return ret_model


