'''
Created on Mar 5, 2014

@author: harrigan
'''

import argparse
import os

import logging as log
from maccelerator import simulate, model
import mdtraj as md
import numpy as np


def run_func(args):

    system, integrator = simulate.deserialize(args.sys_fn, args.int_fn)
    sstate_fn = os.path.join('sstates', 'round-%d.h5' % args.round)
    sstate_traj = md.load(sstate_fn)

    # Pick out the nth frame, loop around
    sstate = sstate_traj[args.traj % sstate_traj.n_frames]

    # Where to save
    traj_out_fn = os.path.join('trajs', 'round-%d' % args.round, 'traj%d.h5' % args.traj)

    # Do it
    simulate.simulate(sstate=sstate, system=system,
                      integrator=integrator, n_spt=args.n_spt,
                      report_stride=args.report, traj_out_fn=traj_out_fn)
    pass

def model_func(args):

    # Load trajectories
    trajs = []
    for round in range(args.round + 1):
        tdir = os.path.join('trajs', 'round-%d' % round)
        trajs += [md.load(os.path.join(tdir, s)) for s in os.listdir(tdir) if s.endswith('.h5')]

    counts, centroids = model.model(trajs=trajs,
                         lagtime=args.lagtime,
                         distance_cutoff=args.distance_cutoff)

    counts_per_state = np.asarray(counts.sum(axis=1)).flatten()
    states_to_sample = np.argsort(counts_per_state)
    if len(states_to_sample) > args.n_tpr:
        states_to_sample = states_to_sample[:args.n_tpr]
    log.info('Generating %d new starting structures.', len(states_to_sample))
    sstates = centroids[states_to_sample]
    sstates.save(os.path.join('sstates', 'round-%d.h5' % (args.round + 1)))



def system_func(args):
    pass



def parse():
    parser = argparse.ArgumentParser(description='Perform accelerated sampling',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sp = parser.add_subparsers()

    parser.add_argument('--sys_fn', help='System xml file', default='../system.xml')
    parser.add_argument('--int_fn', help='Integrator xml file', default='../integrator.xml')

    #===========================================================================
    # Run
    #===========================================================================
    run_p = sp.add_parser('run', help='Run a simulation',
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    run_p.add_argument('--round', '-r', help='Round index', type=int, required=True)
    run_p.add_argument('--traj', '-t', help='Trajectory index', type=int, required=True)
    run_p.add_argument('--n_spt', help='Steps per traj', type=int, required=True)
    run_p.add_argument('--report', help='Report interval', type=int, required=True)

    run_p.set_defaults(func=run_func)


    #===========================================================================
    # Model
    #===========================================================================
    model_p = sp.add_parser('model', help='Create a model and generate new starting structures',
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    model_p.add_argument('--round', '-r', help='Round after which to build model',
                         type=int, required=True)
    model_p.add_argument('--lagtime', '-lt', help='Adaptive modeler lag time',
                         type=int, default=20)
    model_p.add_argument('--distance_cutoff', '-dc', help='Distance cutoff for clustering',
                         type=float, default=0.3)
    model_p.add_argument('--n_tpr', help='Number of new starting states to make',
                         type=int, required=True)
    model_p.set_defaults(func=model_func)


    #===========================================================================
    # System
    #===========================================================================
    system_p = sp.add_parser('system', help='Generate a system of jobs and submission scripts',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    system_p.set_defaults(func=system_func)


    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    log.basicConfig(level=log.INFO)
    parse()
