"""Classes for choosing new starting states."""

import numpy as np
import logging as log
import pickle

__author__ = 'harrigan'


class Adapter(object):
    """Base class for an object that chooses where to start new simulation."""

    def adapt(self, model, params):
        """Return a state from which to start.

        :param params: Simulation parameters.
        """
        raise NotImplementedError

    def seed_states(self, params):
        """Return a state from which to start without requiring a model."""
        raise NotImplementedError

    @property
    def sstatefn(self):
        return "sstate-{round_i}"


class SStates():
    """Objects of this type will be returned from Adapter.adapt."""

    def __init__(self, indices):
        self.indices = indices

    def items(self):
        return self.indices


    def save(self, fn):
        """Save to a file."""
        fn = "{}.pickl".format(fn)
        with open(fn, 'wb') as f:
            pickle.dump(self, f)


class SortCountsAdapter(Adapter):
    """Choose the states from which we've transitioned the fewest times,
    in order
    """


    def adapt(self, model, params):
        """From a counts matrix, pick the best states from which to start.

        :param params: Simulation parameters so we know how many new states
                       to return
        :returns: Indices of states
        """
        counts = model.counts
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

        return SStates(states_to_sample)


class RandomAdapter(Adapter):
    def adapt(self, model, params):
        """Pick random indices

        :param params: So we know how many new states to return
        :returns: Indices of new states
        """
        sstate = SStates(
            np.random.randint(0, model.n_states, params.tpr))

        return sstate


