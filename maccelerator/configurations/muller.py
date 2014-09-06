__author__ = 'harrigan'

import itertools

from simtk import unit
import mdtraj as md
import numpy as np
from simtk import openmm

from msmtoys import muller
from ..simulate import OpenMMSimulator
from ..files import get_fn
from ..model import ClusterModeller
from ..adapt import SortCountsAdapter
from .base import OpenMMConfiguration
from ..convergence import Volume, OpenMMConvergenceChecker
from ..param import AdaptiveParams


def generate_muller_sysint():
    """Set up muller potential."""

    mass = 12.0 * unit.dalton
    temperature = 750 * unit.kelvin
    friction = 100 / unit.picosecond
    timestep = 10.0 * unit.femtosecond

    # Prepare the system
    system = openmm.System()
    mullerforce = muller.MullerForce()
    system.addParticle(mass)
    mullerforce.addParticle(0, [])
    system.addForce(mullerforce)

    # Prepare the integrator
    integrator = openmm.LangevinIntegrator(temperature, friction, timestep)

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


class MullerModeller(ClusterModeller):
    """Model only with x, y"""

    def load_trajs(self, traj_fns):
        """Load only the xy coordinates of an mdtraj trajectory.

        :param traj_fns: List of trajectory filenames.
        """
        trajs = [md.load(tfn).xyz[:, 0, 0:2] for tfn in traj_fns]
        return trajs


class MullerAdapter(SortCountsAdapter):
    # TODO: Generalize
    def adapt(self, params, sstate_out_fn):
        state_indices = super()._adapt(params)
        sstate_positions = self.modeller.clusterer.cluster_centers_[
                           state_indices, :]
        assert sstate_positions.shape[1] == 2
        sstate_positions = np.hstack(
            (sstate_positions, np.zeros((len(sstate_positions), 1))))
        assert sstate_positions.shape[1] == 3

        sstate_traj = make_traj_from_coords(sstate_positions)

        return sstate_traj

    def seed_state(self, params, sstate_out_fn):
        """Start from the bottom right well.

        :param params: Make this many seed states. They will all be the same
        """
        seed_state = make_traj_from_coords([[0.5, 0.0, 0.0]] * params.tpr)
        seed_state.save(sstate_out_fn)
        return seed_state


class MullerParams(AdaptiveParams):
    def __init__(self, spt, tpr, adapt_lt=20, build_lt=20, post_converge=10,
                 run_id=0):
        super().__init__(spt, tpr, adapt_lt, build_lt, post_converge, run_id)


class MullerConfiguration(OpenMMConfiguration):
    def defaults(config):
        """Set defaults."""

        # Simulation report stride
        config.report_stride = 10

        # Whether to perform minimization
        config.minimize = False

        # Start with random initial velocities
        config.random_initial_velocities = True

        # Grid for projection for convergence
        volume = Volume([-1.0, 1.0], [-0.1, 2.0])
        config.grid = Volume.get_grid(volume, resolution=200)

        # Simulator
        config.simulator_class = OpenMMSimulator

        # Modeller
        config.modeller_class = MullerModeller

        # Convergence checker
        config.convchecker_class = OpenMMConvergenceChecker

        # Adapter
        # TODO: Generalize
        config.adapter_class = MullerAdapter

        # Define a function to yield combinations of parameters
        def get_param_grid(run_id):
            spts = [10, 100]
            tprs = [1, 10]

            for spt, tpr in itertools.product(spts, tprs):
                yield MullerParams(spt=spt, tpr=tpr, run_id=run_id)

        config.get_param_grid = get_param_grid


    def __init__(self):
        super().__init__(get_fn('muller_sys.xml'), get_fn('muller_int.xml'))

        self.grid = None
        self.temp = 750
        self.potentialfunc = muller.MullerForce.potential



