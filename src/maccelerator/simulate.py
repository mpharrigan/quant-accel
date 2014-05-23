'''
Created on Mar 5, 2014

@author: harrigan

Code for running rounds of simulation

Code copied from rmcgibbo
'''

import os

from mdtraj.reporters import HDF5Reporter
from simtk.openmm import XmlSerializer, Platform
from simtk.openmm.app import Simulation
from simtk.openmm.app import StateDataReporter

import logging as log
import numpy as np
import mdtraj as md
from mdtraj import io

from quantaccel import tmat_simulation


# My imports
def deserialize(system_xml, integrator_xml):
    """Deserialize system and integrators."""
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

    return system, integrator


def simulate_openmm(sstate, system, integrator, n_spt, report_stride,
                    traj_out_fn,
                    minimize=False, random_initial_velocities=True):
    """Simulate."""

    log.debug('Setting up simulation...')

    # Use cpu platform
    #platform = Platform.getPlatformByName('CPU')
    platform = None

    simulation = Simulation(sstate.topology.to_openmm(), system, integrator,
                            platform)
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
            print "I don't know what temperature to use!!"
            raise

    log.debug('adding reporters...')
    add_reporters(simulation, traj_out_fn, report_stride, n_spt)

    # run dynamics!
    log.debug('Starting dynamics')
    simulation.step(n_spt)

    for reporter in simulation.reporters:
        # explicitly delete the reporters to close file handles
        del reporter


def simulate_muller(args):
    """Load up relevant files and simulate muller using openmm."""

    # Prepare filenames
    sstate_fn, traj_out_fn = get_filenames(args)

    # Load stuff
    system, integrator = deserialize(args.sys_fn, args.int_fn)
    sstate_traj = md.load(sstate_fn)

    # Pick out the nth frame, loop around
    sstate = sstate_traj[args.traj % sstate_traj.n_frames]

    # Do it
    simulate_openmm(sstate=sstate, system=system, integrator=integrator,
                    n_spt=args.n_spt, report_stride=args.report,
                    traj_out_fn=traj_out_fn)


def simulate_tmat(args):
    """Load up relevant files and simulate tmat using quantaccel.TMatSimulator.
    """

    # Prepare filenames
    sstate_fn, traj_out_fn = get_filenames(args)

    # Load stuff
    sim = tmat_simulation.TMatSimulator(args.tmat_fn, get_eigen=False)
    sstate_traj = io.loadh(sstate_fn, 'arr_0')

    # Pick out the nth frame, loop around
    sstate = sstate_traj[args.traj % len(sstate_traj)]

    # Do it
    sim.simulate(state_i=sstate, number_of_steps=args.n_spt, out_fn=traj_out_fn)


def get_filenames(args):
    sstate_fn = os.path.join('sstates', 'round-%d.h5' % args.round)
    traj_out_fn = os.path.join('trajs', 'round-%d' % args.round,
                               'traj%d.h5' % args.traj)
    return sstate_fn, traj_out_fn


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

    def __init__(
            self, reportCallback, reportInterval, total_steps=None, **kwargs):
        super(
            CallbackReporter,
            self).__init__(
            os.devnull,
            reportInterval,
            **kwargs)

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


def add_reporters(simulation, outfn, report_stride, n_spt):
    "Add reporters to a simulation"

    def reporter_callback(report):
        """Callback for processing reporter output"""
        log.debug(report)

    callback_reporter = CallbackReporter(reporter_callback,
                                         report_stride, step=True,
                                         potentialEnergy=True,
                                         temperature=True, time=True,
                                         total_steps=n_spt)

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
    seed = int(
        hashlib.md5(
            '%s%s%s' %
            (plt, os.getpid(), time.time())).hexdigest(), 16)

    return seed % np.iinfo(np.int32).max
