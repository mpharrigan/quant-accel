__author__ = 'harrigan'

import unittest
from unittest import TestCase
import os
from os.path import join as pjoin
import logging

import numpy as np
from numpy.testing import assert_array_equal
from simtk import unit

import maccelerator as maccel
from maccelerator.configurations.muller import (generate_muller_sysint,
                                                make_traj_from_coords)
from maccelerator.simulate import OpenMMSimulator
from maccelerator.testing.test_utils import get_folder



# Disable logging during test
logging.captureWarnings(True)
logging.disable(logging.WARNING)


class TestMullerPrep(TestCase):
    def setUp(self):
        self.tdir = get_folder('mul_prep')

    def test_make(self):
        system, integrator = generate_muller_sysint()

        self.assertFalse(system is None)
        self.assertFalse(integrator is None)
        self.assertEqual(len(system.getForces()), 1)
        self.assertEqual(integrator.getTemperature().value_in_unit(unit.kelvin),
                         750.0)

        sysfn = pjoin(self.tdir, 'muller_sys.xml')
        intfn = pjoin(self.tdir, 'muller_int.xml')

        OpenMMSimulator.serialize(system, integrator, sysfn, intfn)

        self.assertTrue(os.path.exists(pjoin(self.tdir, 'muller_sys.xml')))
        self.assertTrue(os.path.exists(pjoin(self.tdir, 'muller_int.xml')))

        system, integrator = OpenMMSimulator.deserialize(sysfn, intfn)
        self.assertFalse(system is None)
        self.assertFalse(integrator is None)
        self.assertEqual(len(system.getForces()), 1)
        self.assertEqual(integrator.getTemperature().value_in_unit(unit.kelvin),
                         750.0)


    def test_reference(self):
        config = maccel.MullerConfiguration().apply_configuration()
        self.assertEqual(len(config.simulator.system.getForces()), 1)
        self.assertEqual(
            config.simulator.integrator.getTemperature().value_in_unit(
                unit.kelvin), 750.0)


class TestTrajFromNumpy(TestCase):
    def setUp(self):
        self.xy = np.arange(10).reshape(5, 2)
        self.xy = np.hstack((self.xy, np.zeros((5, 1))))

    def test_make_traj_from_coords(self):
        traj = make_traj_from_coords(self.xy)

        self.assertEqual(len(traj.xyz.shape), 3)
        self.assertEqual(traj.xyz.shape[2], 3)
        self.assertEqual(traj.n_frames, 5)
        self.assertEqual(traj.n_atoms, 1)

        assert_array_equal(traj.xyz[:, 0, :], self.xy)


if __name__ == "__main__":
    unittest.main()
