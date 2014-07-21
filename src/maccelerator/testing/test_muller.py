__author__ = 'harrigan'

import unittest
import tempfile
import os
from os.path import join as pjoin

from simtk import unit

import maccelerator as maccel
from maccelerator.configurations.muller import generate_muller_sysint
from maccelerator.simulate import serialize_openmm


class TestMullerPrep(unittest.TestCase):
    def setUp(self):
        self.tdir = tempfile.mkdtemp()

    def test_make(self):
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


if __name__ == "__main__":
    unittest.main()
