'''
Created on Mar 5, 2014

@author: harrigan

Code for running rounds of simulation

Code copied from rmcgibbo
'''

import os
import logging as log

from mdtraj.reporters import HDF5Reporter
from simtk.openmm import XmlSerializer
from simtk.openmm.app import Simulation
from simtk.openmm.app import StateDataReporter
from simtk import openmm
import numpy as np
import mdtraj as md

from mdtraj import io


def generate_openmm_sysint(mass, temperature, friction, timestep):
    """Template for generating openmm files.

    Run this outside of the loop. It is not pickleable so IPython.parallel
    complains
    """

    # Prepare the system
    system = openmm.System()

    # And integrator
    integrator = openmm.LangevinIntegrator(temperature, friction, timestep)

    return system, integrator


def serialize_openmm(system, integrator, sys_fn, int_fn):
    with open(sys_fn, 'w') as sys_f:
        sys_f.write(XmlSerializer.serialize(system))

    with open(int_fn, 'w') as int_f:
        int_f.write(XmlSerializer.serialize(integrator))


class Simulator(object):
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
    def __init__(self, system_xml, integrator_xml, report_stride):
        super(OpenMMSimulator, self).__init__()

        self.report_stride = report_stride
        self.system_xml = system_xml
        self.integrator_xml = integrator_xml
        self.system = None
        self.integrator = None


    def deserialize(self, system_xml, integrator_xml):
        """Deserialize system and integrators.

        :param system_xml: Path to system.xml file
        :param integrator_xml: Path to integrator.xml file
        """

        # load up the system and integrator files
        with open(system_xml) as f:
            system = XmlSerializer.deserialize(f.read())
            # reset the random number seed for any random
            # forces (andersen thermostat, montecarlo barostat)
            for i in range(system.getNumForces()):
                force = system.getForce(i)
                if hasattr(force, 'setRandomNumberSeed'):
                    force.setRandomNumberSeed(random_seed())
        with open(integrator_xml) as f:
            integrator = XmlSerializer.deserialize(f.read())

            # reset the random number seed for a stochastic integrator
            if hasattr(integrator, 'setRandomNumberSeed'):
                integrator.setRandomNumberSeed(random_seed())

        self.system = system
        self.integrator = integrator


    def _simulate(self, sstate, n_steps, traj_out_fn, minimize,
                  random_initial_velocities):
        """Simulate."""

        if self.system is None or self.integrator is None:
            log.debug('Loading system and integrator xmls')
            self.deserialize(self.system_xml, self.integrator_xml)
            log.debug('Xml files loaded. Setting up simulation')
        else:
            log.debug('Xml files were already loaded. Setting up simulation')

        platform = None
        simulation = Simulation(sstate.topology.to_openmm(), self.system,
                                self.integrator, platform)
        simulation.context.setPositions(sstate.openmm_positions(0))
        sanity_check(simulation)

        if minimize:
            log.debug('minimizing...')
            simulation.minimizeEnergy()

        if random_initial_velocities:
            try:
                temp = simulation.integrator.getTemperature()
                simulation.context.setVelocitiesToTemperature(temp)
            except AttributeError:
                raise ValueError("I don't know what temperature to use")

        log.debug('adding reporters...')
        add_reporters(simulation, traj_out_fn, self.report_stride, n_steps)

        # run dynamics!
        log.debug('Starting dynamics')
        simulation.step(n_steps * self.report_stride)

        for reporter in simulation.reporters:
            # explicitly delete the reporters to close file handles
            del reporter

        return True


    @property
    def trajfn(self):
        return "traj-{traj_i}.h5"


#    def simulate_muller(self, args):
#        #TODO: This code should probably be used somewhere
#        """Load up relevant files and simulate muller using openmm."""
#
#        # Prepare filenames
#        sstate_fn, traj_out_fn = get_filenames(args)
#
#        # Load stuff
#        system, integrator = deserialize(args.sys_fn, args.int_fn)
#        sstate_traj = md.load(sstate_fn)
#
#        # Pick out the nth frame, loop around
#        sstate = sstate_traj[args.traj % sstate_traj.n_frames]
#
#        # Do it
#        simulate_openmm(sstate=sstate, system=system, integrator=integrator,
#                        n_spt=args.n_spt, report_stride=args.report,
#                        traj_out_fn=traj_out_fn)


class TMatSimulator(Simulator):
    @property
    def n_states(self):
        """Number of states in the model."""
        return self.t_matrix.shape[0]

    @property
    def trajfn(self):
        return "traj-{traj_i}.h5"

    def __init__(self, t_matrix):
        # Load transition matrix
        #        t_matrix = scipy.io.mmread(tmat_fn)
        #        t_matrix = t_matrix.tocsr()

        self.t_matrix = t_matrix
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

            csr_slicer = slice(
                t_matrix.indptr[sstate],
                t_matrix.indptr[sstate + 1])
            probs = t_matrix.data[csr_slicer]
            colinds = t_matrix.indices[csr_slicer]

            # Find our new state and translate to actual indices
            prob_i = np.sum(np.cumsum(probs) < np.random.rand())
            sstate = colinds[prob_i]

            state_out[i] = sstate

        # Write
        if traj_out_fn is not None:
            io.saveh(traj_out_fn, state_out)
        log.debug('Finished TMat simulation.')
        return state_out


#def get_filenames(args):
#    sstate_fn = os.path.join('sstates', 'round-%d.h5' % args.round)
#    traj_out_fn = os.path.join('trajs', 'round-%d' % args.round,
#                               'traj%d.h5' % args.traj)
#    return sstate_fn, traj_out_fn


def sanity_check(simulation):
    """Lovingly ripped from @rmcgibbo
    """
    positions = simulation.context.getState(
        getPositions=True).getPositions(
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


def add_reporters(simulation, outfn, report_stride, n_spt):
    "Add reporters to a simulation"

    def reporter_callback(report):
        """Callback for processing reporter output"""
        log.debug(report)

    callback_reporter = CallbackReporter(reporter_callback,
                                         report_stride, step=True,
                                         potentialEnergy=True,
                                         temperature=True, time=True,
                                         total_steps=n_spt * report_stride)

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
