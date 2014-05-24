__author__ = 'harrigan'

import os

from simtk import unit
from msmtoys import muller
import mdtraj as md
import numpy as np

from ..simulate import OpenMMSimulator, generate_openmm_sysint
from ..model import ClusterModeller, SortCountsAdapter
from ..configuration import OpenMMConfiguration
from ..param import AdaptiveParams


def generate_muller_sysint():
    """Set up muller potential."""

    mass = 12.0 * unit.dalton
    temperature = 750 * unit.kelvin
    friction = 100 / unit.picosecond
    timestep = 10.0 * unit.femtosecond
    system, integrator = generate_openmm_sysint(mass, temperature,
                                                friction, timestep)

    # Prepare the system
    mullerforce = muller.MullerForce()
    system.addParticle(mass)
    mullerforce.addParticle(0, [])
    system.addForce(mullerforce)

    return system, integrator


class MullerSimulator(OpenMMSimulator):
    def __init__(self, system_xml, integrator_xml):
        super().__init__(system_xml, integrator_xml, report_stride=10)

    def simulate(self, sstate, n_steps, traj_out_fn):
        return super()._simulate(sstate, n_steps, traj_out_fn, minimize=False,
                                 random_initial_velocities=True)


class MullerModeller(ClusterModeller):
    def __init__(self):
        super().__init__()

    def seed_state(self, tpr):
        """Start from the bottom right well.

        :param tpr: Make this many seed states. They will all be the same
        """

        top = md.Topology()
        chain = top.add_chain()
        resi = top.add_residue('NA', chain)
        top.add_atom('C', md.element.carbon, resi)

        xyz = np.array([[[0.5, 0.0, 0.0]]])

        seed_state = md.Trajectory(xyz, top)
        return [seed_state for i in range(tpr)]


class MullerAdapter(SortCountsAdapter):
    pass


class MullerParams(AdaptiveParams):
    @property
    def post_converge(self):
        return 10

    @property
    def adapt_lt(self):
        return 20

    @property
    def build_lt(self):
        return 20


class MullerConfiguration(OpenMMConfiguration):
    def __init__(self, system_xml, integrator_xml):
        super().__init__()

        system_xml = os.path.abspath(system_xml)
        integrator_xml = os.path.abspath(integrator_xml)

        self.simulator = MullerSimulator(system_xml, integrator_xml)
        self.modeller = MullerModeller()
        self.adapter = MullerAdapter(self.modeller)


