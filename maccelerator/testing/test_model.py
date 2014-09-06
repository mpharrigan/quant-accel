__author__ = 'harrigan'
from unittest import TestCase
import unittest
import logging

import numpy as np

import maccelerator as maccel

# Disable logging during test
logging.captureWarnings(True)
logging.disable(logging.WARNING)


class TestTmatModel1(TestCase):
    def make_load_trajs_func(self):
        def load_trajs(self, traj_fns, up_to=None):
            return [np.array([3, 3, 3, 3, 3, 3, 3]),
                    np.array([5, 5, 5, 5, 5, 5])]

        return load_trajs

    def setUp(self):
        config = maccel.AlanineConfiguration()

        maccel.TMatModeller.load_trajs = self.make_load_trajs_func()
        config.apply_configuration()
        self.param = maccel.AlanineParams(0, 0)
        self.model = config.model(None, self.param)
        self.config = config

    def test_convergence(self):
        converge = self.config.check_convergence(self.model, self.param)
        self.assertFalse(converge.converged)

    def test_counts(self):
        should_be = np.zeros((20, 20))
        should_be[3, 3] = 6
        np.testing.assert_array_equal(self.model.counts, should_be)

    def test_tmat(self):
        should_be = np.zeros((20, 20))
        should_be[3, 3] = 1
        np.testing.assert_array_equal(self.model.full_tmat, should_be)

    def test_populations(self):
        should_be = np.zeros(20)
        should_be[3] = 1
        np.testing.assert_array_equal(self.model.full_populations, should_be)


class TestTmatModel2(TestTmatModel1):
    def make_load_trajs_func(self):
        def load_trajs(self, traj_fns, up_to=None):
            return [np.array([3, 3])]
        return load_trajs


if __name__ == "__main__":
    unittest.main()
