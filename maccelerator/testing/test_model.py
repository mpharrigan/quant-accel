__author__ = 'harrigan'
from unittest import TestCase
from os.path import join as pjoin
import os
import unittest
import logging
import pickle

import mdtraj.io
import scipy.io

import maccelerator as maccel
from maccelerator.adapt import SStates
from maccelerator.model import Model
from maccelerator.param import AdaptiveParams
from maccelerator.testing.test_utils import get_folder, get_fn

import numpy as np


# Disable logging during test
logging.captureWarnings(True)
logging.disable(logging.WARNING)


class TestTmatModel(TestCase):
    def setUp(self):
        config = maccel.AlanineConfiguration()

        def load_trajs(self, traj_fns, up_to=None):
            return [np.array([3, 3, 3, 3, 3, 3, 3]),
                    np.array([5, 5, 5, 5, 5, 5])]

        maccel.TMatModeller.load_trajs = load_trajs
        config.apply_configuration()
        self.param = maccel.AlanineParams(0, 0)
        self.model = config.model(None, self.param)
        self.config = config

    def test_convergence(self):
        converge = self.config.check_convergence(self.model, self.param)
        self.assertFalse(converge.converged)


if __name__ == "__main__":
    unittest.main()
