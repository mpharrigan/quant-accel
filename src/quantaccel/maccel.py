#!/usr/bin/env python

"""
Created on Mar 5, 2014

@author: harrigan

Prototype msmaccelerator
"""

import argparse

import logging as log
from maccelerator import simulate, model, generate_system


def run_func(args):
    """Entry point for running a round of simulation."""

    if args.system_type == 'muller':
        simulate.simulate_muller(args)
    elif args.system_type == 'tmat':
        simulate.simulate_tmat(args)
    else:
        raise ValueError('Invalid system type')


def load_trajectories(round_i, load_func):
    raise NotImplementedError('This function has moved')


def load_trajectories_percent(percent, load_func):
    raise NotImplementedError('This function has moved')


def model_func(args):
    """Entry point for modeling and performing adaptive sampling."""

    # Establish differences between system types
    if args.system_type == 'muller':
        model.model_and_adapt_muller(args)
    elif args.system_type == 'tmat':
        model.model_and_adapt_tmat(args)
    else:
        raise ValueError('Invalid system type')

        # TODO: Get rid of distance cutoff as command line?


def dep_system_func(args):
    """Entry point for generating a system of dependent jobs."""
    proj_dir = 'lt-{lagtime}_spt-{n_spt}_tpr-{n_tpr}'.format(**vars(args))
    generate_system.write_dep_jobs(proj_dir, args)


def one_system_init_func(args):
    """Entry point for initializing a one-job run.

    This will set up folder structures and write the bash
    script that actually does the loop.
    """
    proj_dir = 'lt-{lagtime}_spt-{n_spt}_tpr-{n_tpr}'.format(**vars(args))
    generate_system.write_one_job(proj_dir, args)


def one_system_newround_func(args):
    """Entry point for generating one round of jobs."""
    if args.system_type == 'muller':
        args.report = 10
    elif args.system_type == 'tmat':
        args.report = 1
    else:
        raise ValueError("Invalid system type")
    generate_system.write_new_round(args)


def one_system_combi_func(args):
    """Entry point for generating jobs with combinatorical parameters."""

    if args.system_type == 'muller':
        configs = generate_system.create_configs_muller()
    elif args.system_type == 'tmat':
        configs = generate_system.create_configs_tmat()
    else:
        raise ValueError("Invalid system type")

    generate_system.write_combi_jobs(args.runcopy, configs, args.system_type,
                                     args.overwrite)


def parse():
    """Parse arguments."""
    parser = argparse.ArgumentParser(description='Perform accelerated sampling',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sp = parser.add_subparsers()

    parser.add_argument('--system_type', '-st',
                        help="Type of system to run.",
                        choices=['muller', 'tmat'],
                        default='muller')
    parser.add_argument('--sys_fn', help="""System xml file.""",
                        default='../../system.xml')
    parser.add_argument('--int_fn', help='Integrator xml file',
                        default='../../integrator.xml')
    parser.add_argument('--tmat_fn',
                        help='Transition matrix from which to sample',
                        default='../../tProb.mtx')

    #=========================================================================
    # Run
    #=========================================================================
    run_p = sp.add_parser('run', help='Run a simulation',
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    run_p.add_argument('--round', '-r', help='Round index', type=int,
                       required=True)
    run_p.add_argument('--traj', '-t', help='Trajectory index', type=int,
                       required=True)
    run_p.add_argument('--n_spt', help='Steps per traj', type=int,
                       required=True)
    run_p.add_argument('--report', help='Report interval', type=int,
                       required=True)

    run_p.set_defaults(func=run_func)

    #=========================================================================
    # Model
    #=========================================================================
    model_p = sp.add_parser('model',
                            help='Create a model and generate new starting structures',
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    model_p.add_argument('--round', '-r',
                         help='Round after which to build model', type=int,
                         required=True)
    model_p.add_argument('--lagtime', '-lt', help='Adaptive modeler lag time',
                         type=int, default=20)
    model_p.add_argument('--distance_cutoff', '-dc',
                         help='Distance cutoff for clustering', type=float,
                         default=0.2)
    model_p.add_argument('--n_tpr',
                         help='Number of new starting states to make', type=int,
                         required=True)
    model_p.set_defaults(func=model_func)

    #=========================================================================
    # System
    #=========================================================================
    system_p = sp.add_parser('system',
                             help='Generate a system of jobs and submission scripts',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    system_p.add_argument('--n_tpr', help='Number of trajs per round', type=int,
                          required=True)
    system_p.add_argument('--n_spt',
                          help='Number of steps per traj. *Units of rep* ',
                          type=int, required=True)
    system_p.add_argument('--lagtime', '-lt',
                          help='Lagtime at which to make the model.', type=int,
                          required=True)
    system_p.add_argument('--seed_structures',
                          help="""Initial structures to sample from.
                          For openmm-based runs, this should be a trajectory
                          loadable by mdtraj. For tmat-based runs, this should
                          be an HDF5 file containing an array of integers
                          corresponding to seed state indices. It should
                          be loadable with mdtraj.io.loadh()""",
                          default='seed_structures.h5')


    #=========================================================================
    # System -> Dependency
    #=========================================================================
    ssp = system_p.add_subparsers()
    dep_sys_p = ssp.add_parser('dep', help='Use dependency jobs',
                               formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    dep_sys_p.add_argument('--n_tpj', help='Number of trajectories per job',
                           type=int, default=1)
    dep_sys_p.add_argument('--n_round', help='Number of rounds', type=int,
                           default=20)
    dep_sys_p.add_argument('--start_from', '-sf', help="""Start at a particular round
                          Note: ! This doesn't do any fancy error checking.""",
                           type=int)
    dep_sys_p.set_defaults(func=dep_system_func)

    #=========================================================================
    # System -> One job
    #=========================================================================
    one_sys_p = ssp.add_parser('one', help='Put everything in one job',
                               formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    #===========================================================================
    # System -> One job -> Init
    #===========================================================================
    osp = one_sys_p.add_subparsers()
    one_init_p = osp.add_parser('init', help="Initialize directory",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    one_init_p.add_argument('--n_postconverge',
                            help="Number of rounds to run post convergence",
                            default=10, type=int)
    one_init_p.add_argument('--how', help="'how' string for building msm et al",
                            default='rnew')
    one_init_p.add_argument('--version',
                            help="'version' number mk{version} for building msm et al.",
                            default=7, type=int)
    one_init_p.set_defaults(func=one_system_init_func)

    #===========================================================================
    # System -> One job -> New round
    #===========================================================================
    one_newround_p = osp.add_parser('newround',
                                    help="Write jobs for a new round",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    one_newround_p.add_argument('--round_i', help="Round to build stuff for",
                                type=int, required=True)
    one_newround_p.set_defaults(func=one_system_newround_func)

    #===========================================================================
    # Combi
    #===========================================================================
    one_combi = sp.add_parser('combi', help="Make combinatorical jobs",
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    one_combi.add_argument('--runcopy', help="Index of the run to generate",
                           required=True, type=int)
    one_combi.add_argument('--overwrite', help="Overwrite previous runcopy",
                           action='store_true')
    one_combi.set_defaults(func=one_system_combi_func, overwrite=False)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    log.basicConfig(level=log.INFO)
    parse()
