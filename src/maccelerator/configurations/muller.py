__author__ = 'harrigan'

import os

from simtk import unit
from msmtoys import muller
import mdtraj as md
import numpy as np

from ..simulate import OpenMMSimulator, generate_openmm_sysint
from ..model import ClusterModeller, SortCountsAdapter
from ..configuration import OpenMMConfiguration
from ..check_convergence import PopulationProjectionTVD, Volume
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


def make_traj_from_coords(xyz):
    """Take numpy array and turn it into a one particle trajectory

    :param xyz: (n_frames, 3) np.ndarray
    """

    top = md.Topology()
    chain = top.add_chain()
    resi = top.add_residue('NA', chain)
    top.add_atom('C', md.element.carbon, resi)

    xyz = np.asarray(xyz)
    xyz = xyz[:, np.newaxis, :]
    traj = md.Trajectory(xyz, top)
    return traj


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
        seed_state = make_traj_from_coords([[0.5, 0.0, 0.0]] * params.tpr)
        seed_state.save('{}.h5'.format(params.seed_state_fn))
        return seed_state

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


class MullerConvchecker(PopulationProjectionTVD):
    def __init__(self, modeller):
        volume = Volume([-1.5, 1.2], [-0.2, 3.0])
        grid = Volume.get_grid(volume, resolution=200)

        super().__init__(modeller, threshold=0.6, grid=grid,
                         potentialfunc=muller.MullerForce.potential, temp=750)


class MullerAdapter(SortCountsAdapter):
    def adapt(self, params):
        state_indices = super().adapt(params)
        sstate_positions = self.modeller.clusterer.cluster_centers_[
                           state_indices, :]
        assert sstate_positions.shape[1] == 2
        sstate_positions = np.hstack((sstate_positions,
                                      np.zeros((len(sstate_positions), 1))))
        assert sstate_positions.shape[1] == 3

        sstate_traj = make_traj_from_coords(sstate_positions)
        #TODO Save
        return sstate_traj


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
        self.convchecker = MullerConvchecker(self.modeller)
        self.adapter = MullerAdapter(self.modeller)


