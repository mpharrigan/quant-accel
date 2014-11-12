"""
Created on Mar 5, 2014

@author: harrigan

Code for running rounds of simulation

Code copied from rmcgibbo
"""

import os
import logging

from mdtraj.reporters import HDF5Reporter
from simtk.openmm import XmlSerializer
from simtk.openmm.app import Simulation, StateDataReporter
import scipy.sparse
import numpy as np

from mdtraj import io


log = logging.getLogger(__name__)


class Simulator:
    def __init__(self, config):
        pass

    def simulate(self, sstate, params, traj_out_fn):
        """Run a simulation

        :param sstate: Starting state
        :param params: Parameter object that contains number of steps to take
        :param traj_out_fn: Where to save the trajectory
        """
        raise NotImplementedError

    @property
    def trajfn(self):
        return "traj-{traj_i}"


class OpenMMSimulator(Simulator):
    def __init__(self, config):
        super().__init__(config)

        self.report_stride = config.report_stride
        self.system = config.system
        self.integrator = config.integrator

        self.minimize = config.minimize
        self.random_initial_velocities = config.random_initial_velocities

    @classmethod
    def serialize(cls, system, integrator, sys_fn, int_fn):
        """Serialize Openmm system and integrator to files."""
        with open(sys_fn, 'w') as sys_f:
            sys_f.write(XmlSerializer.serialize(system))

        with open(int_fn, 'w') as int_f:
            int_f.write(XmlSerializer.serialize(integrator))

    @classmethod
    def deserialize(cls, system_xml, integrator_xml):
        """Deserialize system and integrators.

        :param system_xml: Path to system.xml file
        :param integrator_xml: Path to integrator.xml file
        """

        # Load system
        with open(system_xml, 'r') as f:
            system = XmlSerializer.deserialize(f.read())

        # reset the random number seed for any random forces
        for i in range(system.getNumForces()):
            force = system.getForce(i)
            try:
                force.setRandomNumberSeed(random_seed())
            except AttributeError:
                pass

        # Load integrator
        with open(integrator_xml, 'r') as f:
            integrator = XmlSerializer.deserialize(f.read())

        # reset the random number seed for a stochastic integrator
        try:
            integrator.setRandomNumberSeed(random_seed())
        except AttributeError:
            pass

        return system, integrator

    def simulate(self, sstate, params, traj_out_fn):
        """Run a simulation

        :param sstate: Starting state
        :param params: Parameter object that contains number of steps to take
        :param traj_out_fn: Where to save the trajectory
        """

        simulation = Simulation(sstate.topology.to_openmm(),
                                self.system, self.integrator)
        simulation.context.setPositions(sstate.openmm_positions(0))
        sanity_check(simulation)

        if self.minimize:
            log.debug('minimizing...')
            simulation.minimizeEnergy()

        if self.random_initial_velocities:
            try:
                temp = simulation.integrator.getTemperature()
                simulation.context.setVelocitiesToTemperature(temp)
            except AttributeError:
                raise ValueError("I don't know what temperature to use")

        log.debug('adding reporters...')
        add_reporters(simulation, traj_out_fn, self.report_stride,
                      params.spt * self.report_stride)

        # Run dynamics!
        log.debug('Starting dynamics')
        simulation.step(params.spt * self.report_stride)

        for reporter in simulation.reporters:
            # explicitly delete the reporters to close file handles
            del reporter

        return True

    @property
    def trajfn(self):
        return "traj-{traj_i}.h5"


class TMatSimulator(Simulator):
    @property
    def n_states(self):
        """Number of states in the model."""
        return self.t_matrix.shape[0]

    @property
    def trajfn(self):
        return "traj-{traj_i}.h5"

    def __init__(self, config):
        super().__init__(config)

        # TODO: Change to sampling from dense matrix
        self.t_matrix = scipy.sparse.csr_matrix(config.ref_msm.transmat_)
        log.info('Using transition matrix of shape %s', self.t_matrix.shape)

    def simulate(self, sstate, params, traj_out_fn):
        """We run some KMC dynamics, and then send back the results.

        :param sstate: Initial state index
        :param params: Contains length of trajectory to return. Note:
                        We actually take one fewer *step* because we include
                        the initial state in our trajectory
        :param traj_out_fn: Optionally write out a trajectory in HDF5 format
        """
        log.debug('Starting TMat simulation...')

        n_steps = params.spt
        t_matrix = self.t_matrix
        state_out = np.zeros(n_steps, dtype=int)

        # Set first state to initial state
        state_out[0] = sstate

        for i in range(1, n_steps):
            # Get stuff from our sparse matrix

            csr_slicer = slice(t_matrix.indptr[sstate],
                               t_matrix.indptr[sstate + 1])
            probs = t_matrix.data[csr_slicer]
            colinds = t_matrix.indices[csr_slicer]

            # Find our new state and translate to actual indices
            prob_i = np.sum(np.cumsum(probs) < np.random.rand())
            sstate = colinds[prob_i]

            state_out[i] = sstate

        # Write
        if traj_out_fn is not None:
            io.saveh(traj_out_fn, state_traj=state_out)
        log.debug('Finished TMat simulation.')
        return state_out


def sanity_check(simulation):
    """Lovingly ripped from @rmcgibbo
    """
    positions = simulation.context.getState(getPositions=True).getPositions(
        asNumpy=True)
    for atom1, atom2 in simulation.topology.bonds():
        d = np.linalg.norm(
            positions[atom1.index, :] - positions[atom2.index, :])
        if not d < 0.3:
            log.error(positions[atom1.index, :])
            log.error(positions[atom2.index, :])
            raise ValueError(
                'atoms are bonded according to topology but not close by '
                'in space: %s. %s' % (d, positions))


class CallbackReporter(StateDataReporter):
    """An openmmm reporter subclass to send report to a callback function."""

    def __init__(self, reportCallback, reportInterval, total_steps=None,
                 **kwargs):

        self.f = open(os.devnull, 'wb')
        super().__init__(self.f, reportInterval, **kwargs)

        self.total_steps = total_steps
        self.reportCallback = reportCallback
        self.headers = None

    def report(self, simulation, state):
        if not self._hasInitialized:
            self._initializeConstants(simulation)
            self.headers = self._constructHeaders()
            self._hasInitialized = True

        # Check for errors.
        self._checkForErrors(simulation, state)

        # Query for the values
        values = self._constructReportValues(simulation, state)

        content = dict(zip(self.headers, values))

        if self.total_steps is not None and 'Step' in content:
            progress = (100 * content['Step'] / self.total_steps)
            content['Progress (%s)'] = '%.1f%%' % progress

        self.reportCallback(content)

    def __del__(self):
        self.f.close()


def add_reporters(simulation, outfn, report_stride, total_steps):
    """Add reporters to a simulation"""

    def reporter_callback(report):
        """Callback for processing reporter output"""
        log.debug(report)

    callback_reporter = CallbackReporter(reporter_callback, report_stride,
                                         step=True, potentialEnergy=True,
                                         temperature=True, time=True,
                                         total_steps=total_steps)

    h5_reporter = HDF5Reporter(outfn, report_stride, coordinates=True,
                               time=True, cell=True, potentialEnergy=True,
                               kineticEnergy=True, temperature=True)

    simulation.reporters.append(callback_reporter)
    simulation.reporters.append(h5_reporter)


def random_seed():
    """Get a seed for a random number generator, based on the current platform,
    pid, and wall clock time.

    Returns
    -------
    seed : int
        The seed is a 32-bit int
    """
    import platform
    import time
    import hashlib

    plt = ''.join(platform.uname())
    seed_str = '%s%s%s' % (plt, os.getpid(), time.time())
    seed_str = seed_str.encode('utf-8')
    seed = int(hashlib.md5(seed_str).hexdigest(), 16)

    return seed % np.iinfo(np.int32).max
