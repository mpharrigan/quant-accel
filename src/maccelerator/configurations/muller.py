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

    def simulate(self, sstate, params, traj_out_fn):
        return super()._simulate(sstate, params.spt, traj_out_fn,
                                 minimize=False,
                                 random_initial_velocities=True)


class MullerModeller(ClusterModeller):
    def __init__(self):
        super().__init__()

    def seed_state(self, params):
        """Start from the bottom right well.

        :param params: Make this many seed states. They will all be the same
        """

        top = md.Topology()
        chain = top.add_chain()
        resi = top.add_residue('NA', chain)
        top.add_atom('C', md.element.carbon, resi)

        xyz = np.array([[[0.5, 0.0, 0.0]]])

        seed_state = md.Trajectory(xyz, top)
        return [seed_state for _ in range(params.tpr)]

    def model(self, traj_fns, params):
        trajs = MullerModeller.load_xy(traj_fns)
        super()._model(trajs, lagtime=params.adapt_lt)

    @staticmethod
    def load_xy(traj_fns):
        """Load only the xy coordinates of an mdtraj trajectory.

        :param traj_fns: List of trajectory filenames.
        """
        trajs = [md.load(tfn).xyz[:, 0, 0:2] for tfn in traj_fns]
        return trajs


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


