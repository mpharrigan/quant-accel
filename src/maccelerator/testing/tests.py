__author__ = 'harrigan'

import unittest
import tempfile
import pickle
import os
from os.path import join as pjoin

import numpy as np
import scipy.sparse
import mdtraj as md

import maccelerator as maccel
from maccelerator.testing.utils import get_fn
from simtk import unit


class TestMullerPrep(unittest.TestCase):
    def setUp(self):
        self.tdir = tempfile.mkdtemp()

    def test_make(self):
        from maccelerator.configurations.muller import generate_muller_sysint
        from maccelerator.simulate import serialize_openmm

        system, integrator = generate_muller_sysint()

        self.assertFalse(system is None)
        self.assertFalse(integrator is None)
        self.assertEqual(len(system.getForces()), 1)
        self.assertEqual(integrator.getTemperature().value_in_unit(unit.kelvin),
                         750.0)

        sysfn = pjoin(self.tdir, 'muller_sys.xml')
        intfn = pjoin(self.tdir, 'muller_int.xml')

        serialize_openmm(system, integrator, sysfn, intfn)

        self.assertTrue(os.path.exists(pjoin(self.tdir, 'muller_sys.xml')))

        config = maccel.MullerConfiguration(sysfn, intfn)
        config.simulator.deserialize(config.simulator.system_xml,
                                     config.simulator.integrator_xml)
        self.assertEqual(len(config.simulator.system.getForces()), 1)
        self.assertEqual(
            config.simulator.integrator.getTemperature().value_in_unit(
                unit.kelvin), 750.0)


class TestMullerSampling(unittest.TestCase):
    def setUp(self):
        self.config = maccel.MullerConfiguration(
            get_fn('muller_sys.xml'), get_fn('muller_int.xml'))
        self.tmpdir = tempfile.mkdtemp()

    def test_sstate(self):
        params = maccel.configurations.MullerParams(spt=100, tpr=20, run_id=-1)
        params.rundir = self.tmpdir
        params.make_directories()
        params._round_i = 5

        sstate = self.config.modeller.seed_state(params)
        for ss in sstate:
            np.testing.assert_array_equal(ss.xyz, [[[0.5, 0.0, 0.0]]])

        sstate_load = md.load(pjoin(self.tmpdir, 'sstates', 'seed-6.h5'))
        for ss in sstate_load:
            np.testing.assert_array_equal(ss.xyz, [[[0.5, 0.0, 0.0]]])

    def test_sampling_length(self):
        params = maccel.configurations.MullerParams(spt=100, tpr=1, run_id=-1)
        sstate = self.config.modeller.seed_state(params)

        success = self.config.simulator.simulate(sstate[0], params,
                                                 pjoin(self.tmpdir, 'trj1.h5'))

        self.assertTrue(success)

        traj = md.load(pjoin(self.tmpdir, 'trj1.h5'))
        self.assertEqual(traj.n_frames, 100)

        # Test vectorization
        traj_xy = maccel.configurations.muller.MullerModeller.load_xy(
            [pjoin(self.tmpdir, 'trj1.h5')])
        self.assertEqual(len(traj_xy), 1)
        s1, s2 = traj_xy[0].shape
        self.assertEqual(s1, 100)
        self.assertEqual(s2, 2)


class TestMullerRun(unittest.TestCase):
    def setUp(self):
        configuration = maccel.MullerConfiguration(
            get_fn('muller_sys.xml'), get_fn('muller_int.xml'))
        param = maccel.SimpleParams(spt=5000, tpr=2)
        rundir = pjoin(tempfile.mkdtemp())
        self.rundir = rundir
        print(self.rundir)
        self.run = maccel.MAccelRun(configuration, param, rundir)

    def test_run(self):
        pass
        self.run.run()

        #TODO: Add asserts


class TestAlanineRun(unittest.TestCase):
    def setUp(self):
        configuration = maccel.AlanineConfiguration(get_fn('ala.msm.pickl'),
                                                    get_fn('ala.centers.h5'))
        param = maccel.SimpleParams(spt=50, tpr=2)
        rundir = pjoin(tempfile.mkdtemp(), 'runz')
        self.rundir = rundir
        self.run = maccel.MAccelRun(configuration, param, rundir)

    def test_run(self):
        self.run.run()
        #TODO: Add asserts


class TestTMatSampling(unittest.TestCase):
    def setUp(self):
        tmat = scipy.sparse.csr_matrix(
            np.diag(np.ones(5), 1)  # 6 state linear transition model
        )
        self.simulator = maccel.simulate.TMatSimulator(tmat)
        self.tmpdir = tempfile.mkdtemp()

    def test_fn(self):
        self.assertEqual(self.simulator.trajfn.format(traj_i=59), 'traj-59.h5')

    def test_length(self):
        params = maccel.param.AdaptiveParams(spt=6, tpr=1)
        seq = self.simulator.simulate(sstate=0, params=params,
                                      traj_out_fn=pjoin(self.tmpdir, 'trj1.h5'))

        self.assertEqual(len(seq), 6)
        np.testing.assert_array_equal(seq, range(6))

        self.assertTrue(os.path.exists(pjoin(self.tmpdir, 'trj1.h5')))



