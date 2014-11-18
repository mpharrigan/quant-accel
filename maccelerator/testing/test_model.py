__author__ = 'harrigan'
from unittest import TestCase
import unittest
import logging

import numpy as np
import maccelerator as maccel


logging.captureWarnings(True)
logging.disable(logging.WARNING)


class TestTmatModel_LongTransitions(TestCase):
    """Test model building with a trajectory that visits all states

    0 --> 19
    19 --> 0
    """

    def make_load_trajs_func(self):
        def load_trajs(self, traj_fns, up_to=None):
            return [np.arange(20),
                    19 - np.arange(20)]

        return load_trajs

    def setUp(self):
        config = maccel.AlanineConfiguration()

        maccel.TMatModeller.load_trajs = self.make_load_trajs_func()
        config.apply_configuration()
        self.param = maccel.AlanineParams(0, 0)
        self.model = config.model(None, self.param)
        self.config = config

    def test_convergence(self):
        # None of these tests should be converged
        # They should not throw errors, however
        converge = self.config.check_convergence(self.model, self.param)
        self.assertFalse(converge.converged)

    def test_counts(self):
        # Expect one diagonal on the top and bottom
        should_be = np.diag(np.ones(19), 1) + np.diag(np.ones(19), -1)

        np.testing.assert_array_equal(self.model.full_counts, should_be)
        np.testing.assert_array_equal(self.model.adapt_counts, should_be)
        np.testing.assert_array_equal(self.model.found_states, np.arange(20))

    def test_tmat(self):
        # Normalized version of the counts. Special care is taken
        # for the first and last state
        should_be = (np.diag(np.ones(19), 1) + np.diag(np.ones(19), -1)) * 0.5
        should_be[0, 1] = 1
        should_be[-1, -2] = 1
        np.testing.assert_array_equal(self.model.tmat, should_be)

    def test_populations(self):
        # As per above, the first and last state are special
        should_be = np.ones(20) * 2
        should_be[0] = 1
        should_be[-1] = 1
        should_be /= np.sum(should_be)
        np.testing.assert_array_equal(self.model.populations, should_be)


class TestTmatModel_Transitions(TestTmatModel_LongTransitions):
    """Test with two short trajectories."""

    def make_load_trajs_func(self):
        def load_trajs(self, traj_fns, up_to=None):
            return [np.array([2, 3, 3, 2, 2]),
                    np.array([3, 2, 2, 3, 3])]

        # Counts:
        # 2  2
        # 2  2

        return load_trajs

    def test_counts(self):
        should_be_adapt = np.array([[2, 2], [2, 2]])
        should_be_full = np.zeros((20, 20))
        should_be_full[2:4, 2:4] = should_be_adapt

        np.testing.assert_array_equal(self.model.full_counts, should_be_full)
        np.testing.assert_array_equal(self.model.adapt_counts, should_be_adapt)
        np.testing.assert_array_equal(self.model.found_states, np.array([2, 3]))

    def test_tmat(self):
        should_be = np.zeros((20, 20))
        should_be[2:4, 2:4] = 0.5 * np.ones((2, 2))
        np.testing.assert_array_equal(self.model.tmat, should_be)

    def test_populations(self):
        should_be = np.zeros(20)
        should_be[2:4] = 0.5
        np.testing.assert_array_equal(self.model.populations, should_be)


class TestTmatModel_LongStay(TestTmatModel_LongTransitions):
    """Test ergodic trimming and 1 state models by setting up two
    trajectories that just stay in their state.
    """

    def make_load_trajs_func(self):
        def load_trajs(self, traj_fns, up_to=None):
            return [np.array([3, 3, 3, 3, 3, 3, 3]),
                    np.array([5, 5, 5, 5, 5, 5])]

        return load_trajs


    def test_counts(self):
        adapt_counts = np.array([
            [6, 0],
            [0, 5]
        ])
        full_counts = np.zeros((20, 20))
        full_counts[3, 3] = 6
        full_counts[5, 5] = 5
        np.testing.assert_array_equal(self.model.full_counts, full_counts)
        np.testing.assert_array_equal(self.model.adapt_counts, adapt_counts)
        np.testing.assert_array_equal(self.model.found_states, np.array([3, 5]))

    def test_tmat(self):
        should_be = np.zeros((20, 20))
        should_be[3, 3] = 1
        np.testing.assert_array_equal(self.model.tmat, should_be)

    def test_populations(self):
        should_be = np.zeros(20)
        should_be[3] = 1
        np.testing.assert_array_equal(self.model.populations, should_be)


class TestTmatModel_ShortStay(TestTmatModel_LongStay):
    """Also test the shortest possible trajectory."""

    def make_load_trajs_func(self):
        def load_trajs(self, traj_fns, up_to=None):
            return [np.array([3, 3])]

        return load_trajs

    def test_counts(self):
        should_be = np.zeros((20, 20))
        should_be[3, 3] = 1
        np.testing.assert_array_equal(self.model.full_counts, should_be)
        np.testing.assert_array_equal(self.model.adapt_counts, [[1]])
        np.testing.assert_array_equal(self.model.found_states, np.array([3]))


if __name__ == "__main__":
    unittest.main()
