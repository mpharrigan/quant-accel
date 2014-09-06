"""Classes for choosing new starting states."""

import logging
import pickle

import numpy as np


__author__ = 'harrigan'

log = logging.getLogger(__name__)


class Adapter:
    """Base class for an object that chooses where to start new simulation."""

    def __init__(self, config):
        # pylint: disable=unused-argument
        pass

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

    @classmethod
    def load(cls, fn):
        with open(fn, 'rb') as f:
            return pickle.load(f)


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
        return SStates(np.random.randint(0, model.tot_n_states, params.tpr))

    def seed_states(self, params):
        return SStates([0] * params.tpr)


