"""Configure types of runs."""

from simtk import openmm, unit
import logging as log


class Configuration(object):
    def __init__(self):
        self.simulator = None
        self.modeler = None

    @staticmethod
    def get_param_grid():
        raise NotImplementedError


class OpenMMConfiguration(Configuration):
    def generate(self, mass, temperature, friction, timestep):
        """Template for generating openmm files."""

        # Prepare the system
        system = openmm.System()

        # And integrator
        integrator = openmm.LangevinIntegrator(temperature, friction, timestep)

        return system, integrator

    def serialize(self, system, integrator, out_sys_fn, out_int_fn):
        #TODO: This code should probably go elsewhere

        with open(out_sys_fn, 'w') as f:
            f.write(openmm.XmlSerializer.serialize(system))
        with open(out_int_fn, 'w') as f:
            f.write(openmm.XmlSerializer.serialize(integrator))


class TMatConfiguration(Configuration):
    pass





