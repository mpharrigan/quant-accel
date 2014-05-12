__author__ = 'harrigan'

from maccelerator.configuration import base_classes as bc
import maccelerator as maccel
from simtk import unit
from msmtoys import muller
import mdtraj as md
import numpy as np


class MullerSimulator(maccel.OpenMMSimulator):
    def __init__(self):
        super().__init__(report_stride=10)

    def simulate(self, sstate, n_steps):
        super().simulate(sstate, n_steps, minimize=False,
                         random_initial_velocities=True)


class MullerModeller(maccel.ClusterModeller):
    def __init__(self):
        super().__init__()

    def seed_state(self):
        """Start from the bottom right well."""

        top = md.Topology()
        chain = top.add_chain()
        resi = top.add_residue('NA', chain)
        top.add_atom('C', md.element.carbon, resi)

        xyz = np.array([[[0.5, 0.0, 0.0]]])

        seed_state = md.Trajectory(xyz, top)
        return seed_state


class MullerConfiguration(bc.OpenMMConfiguration):
    def __init__(self):
        super().__init__()

        self.simulator = MullerSimulator()
        self.modeler = MullerModeller()


    def generate(self):
        """Set up muller potential."""

        mass = 12.0 * unit.dalton
        temperature = 750 * unit.kelvin
        friction = 100 / unit.picosecond
        timestep = 10.0 * unit.femtosecond
        system, integrator = super().generate(mass, temperature, friction,
                                              timestep)

        # Prepare the system
        mullerforce = muller.MullerForce()
        system.addParticle(mass)
        mullerforce.addParticle(0, [])
        system.addForce(mullerforce)

        return system, integrator
