#!/usr/bin/env python

'''
Created on Mar 5, 2014

@author: harrigan
'''

import argparse
import os
import stat
import shutil

import logging as log
from maccelerator import simulate, model
import mdtraj as md
import numpy as np

PBS_HEADER = """
#PBS -l nodes=1:ppn=1
#PBS -l walltime={hours}:00:00
#PBS -l mem=8gb
#PBS -j oe
#PBS -M harrigan@stanford.edu
#PBS -m ae


cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1

"""

SIMULATE_JOB = PBS_HEADER + """
for i in {{{start_i}..{end_i}}}
do
    maccel.py run --round {round_i} --traj $i --n_spt {n_spt} --report {report} &> jobs/{job_fn}.log
done
"""

MODEL_JOB = PBS_HEADER + """
maccel.py model --round {round_i} --lagtime {lagtime} --n_tpr {n_tpr} &> jobs/{job_fn}.log
"""

SIMULATE_SUBMIT = """S{traj_i}=`qsub {dep} jobs/{job_fn}.job`
if echo "$S{traj_i}" | grep -qi "invalid credential"; then echo "error: {job_fn}; exit 1; fi"""

MODEL_SUBMIT = """M{round_i}=`qsub {dep} jobs/{job_fn}.job`
if echo "$M{traj_i}" | grep -qi "invalid credential"; then echo "error: {job_fn}; exit 1; fi"""


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

def load_trajectories(round_i):
    # Load trajectories
    trajs = []
    for cround in range(round_i + 1):
        tdir = os.path.join('trajs', 'round-%d' % cround)
        trajs += [md.load(os.path.join(tdir, s)) for s in os.listdir(tdir) if s.endswith('.h5')]
        
    # Stats
    traj_len = trajs[0].n_frames
    wall_steps = traj_len * (round_i + 1)
    return wall_steps, trajs

def model_func(args):

    _, trajs = load_trajectories(args.round)
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
    proj_dir = 'lt-{lagtime}_spt-{n_spt}_tpr-{n_tpr}'.format(**vars(args))
    start_from = args.start_from
    
    make_dir = start_from is None
    last_id = None
    
    if make_dir:
        # Make project
        log.info('Creating directory %s', proj_dir)
        os.mkdir(proj_dir)
    
        # Make starting structures
        os.mkdir(os.path.join(proj_dir, 'sstates'))
        shutil.copy(args.seed_structures, os.path.join(proj_dir, 'sstates', 'round-0.h5'))
    
        # Make jobs dir
        os.mkdir(os.path.join(proj_dir, 'jobs'))
        
        os.mkdir(os.path.join(proj_dir, 'trajs'))
        
        start_from = 0
    else:
        try:
            with open(os.path.join(proj_dir, 'last.id'), 'r') as f:
                last_id = f.readline().strip()
        except IOError:
            pass
    
    
    submit_lines = []
    for round_i in range(start_from, args.n_round):
        os.mkdir(os.path.join(proj_dir, 'trajs', 'round-%d' % round_i))

        m_dep = '-W depend=afterok'

        # Make simulate jobs
        for traj_i in range(0, args.n_tpr, args.n_tpj):
            job_fn = 'round-{round_i}_traj-{traj_i}'.format(round_i=round_i,
                                                                traj_i=traj_i)
            
            with open(os.path.join(proj_dir, 'jobs', "%s.job" % job_fn), 'w') as job_f:
                job_f.write(SIMULATE_JOB.format(round_i=round_i,
                                                n_spt=args.n_spt * args.report,
                                                report=args.report, job_fn=job_fn,
                                                hours=1,
                                                start_i = traj_i,
                                                end_i = traj_i + args.n_tpj - 1))
                if round_i > start_from:
                    dep = '-W depend=afterok:$M{pti}'.format(pti=round_i - 1)
                elif round_i == start_from and last_id is not None:
                    dep = '-W depend=afterok:{last_id}'.format(last_id=last_id)
                else:
                    dep = ''

                submit_lines += [SIMULATE_SUBMIT.format(traj_i=traj_i,
                                                      job_fn=job_fn,
                                                      dep=dep)]
                m_dep += ':$S{traj_i}'.format(traj_i=traj_i)

        # Make model job
        job_fn = 'round-{round_i}_model'.format(round_i=round_i)
        with open(os.path.join(proj_dir, 'jobs', "%s.job" % job_fn), 'w') as job_f:
            job_f.write(MODEL_JOB.format(round_i=round_i,
                                         lagtime=args.lagtime,
                                         n_tpr=args.n_tpr, hours=3,
                                         job_fn=job_fn))
            submit_lines += [MODEL_SUBMIT.format(job_fn=job_fn,
                                                dep=m_dep,
                                                round_i=round_i)]


    submit_lines += ['echo $M{maxround} > last.id'.format(maxround=args.n_round-1)]
    with open(os.path.join(proj_dir, 'submit.sh'), 'w') as sub_f:
        sub_f.write('\n'.join(submit_lines))

    st = os.stat(os.path.join(proj_dir, 'submit.sh'))
    os.chmod(os.path.join(proj_dir, 'submit.sh'), st.st_mode | stat.S_IEXEC)



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
                         type=float, default=0.2)
    model_p.add_argument('--n_tpr', help='Number of new starting states to make',
                         type=int, required=True)
    model_p.set_defaults(func=model_func)


    #===========================================================================
    # System
    #===========================================================================
    system_p = sp.add_parser('system', help='Generate a system of jobs and submission scripts',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    system_p.add_argument('--n_tpr', help='Number of trajs per round', type=int,
                          required=True)
    system_p.add_argument('--n_spt', help='Number of steps per traj. *Units of rep* ',
                          type=int, required=True)
    system_p.add_argument('--lagtime', '-lt', help='Lagtime at which to make the model.',
                          type=int, required=True)
    system_p.add_argument('--n_round', help='Number of rounds',
                          type=int, default=20)
    system_p.add_argument('--seed_structures', help='Initial structures',
                          default='seed_structures.h5')
    system_p.add_argument('--report', help='Report interval', type=int,
                          default=10)
    system_p.add_argument('--start_from', '-sf', help="""Start at a particular round
                          Note: ! This doesn't do any fancy error checking.""",
                          type=int)
    system_p.add_argument('--n_tpj', help='Number of trajectories per job',
                          type=int, default=1)
    system_p.set_defaults(func=system_func)


    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    log.basicConfig(level=log.INFO)
    parse()
