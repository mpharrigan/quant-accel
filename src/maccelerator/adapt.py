"""Classes for choosing new starting states."""

import numpy as np
import logging as log

__author__ = 'harrigan'


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
        :returns: Indices of states
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


class RandomAdapter(Adapter):
    def __init__(self, modeller):
        super().__init__()
        self.modeller = modeller

    def adapt(self, params):
        """Pick random indices

        :param params: So we know how many new states to return
        :returns: Indices of new states
        """
        return np.random.randint(0, self.modeller.msm.n_states, params.tpr)


