"""
Get implied timescales from the gold trajectory so we can get
    - a sense of how big a step is (lag time to use for other models)
    - gold timescales

Trajs + Cluster -> Assignments

Do a bunch of times:
    Assignments + Lag time -> MSM
    MSM -> Implied timescales

Is it worth using stock MSMBuilder stuff? Probably not, esp. to do
kmeans.

"""

import matplotlib
matplotlib.use('Agg')

import argparse
import logging as log
from msmbuilder import clustering, MSMLib as msmlib, msm_analysis as msma
from quantaccel import toy
import os
import mdtraj
import numpy as np
import stat
from matplotlib import pyplot as pp
import re

CLUSTER_JOB_TEMPLATE = """
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -l mem=8gb
#PBS -j oe
#PBS -o .
#PBS -M harrigan@stanford.edu
#PBS -m bea

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8

python {src_dir}/quantaccel/get_implied_timescales.py cluster &> cluster.log
"""

IT_JOB_TEMPLATE = """
#PBS -l nodes=1:ppn=1
#PBS -l walltime=6:00:00
#PBS -l mem=8gb
#PBS -j oe
#PBS -o .
#PBS -M harrigan@stanford.edu
#PBS -m a

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1

python {src_dir}/quantaccel/get_implied_timescales.py it {lagtime} &> lagtime-{lagtime}.log
"""

SUBMIT_TEMPLATE = \
"""#!/bin/bash

{qsubs}

"""

QSUB_TEMPLATE = \
"""qsub lagtime-{lagtime}.job
sleep 0.5
"""

def cluster_func(args):
    """Cluster trajectories

    Properties of args:
        :traj_dir: Load trajectories from args.traj_dir
        :load_stride: Stride by args.load_stride.
        :distance_cutoff: Cutoff for number of clusters
        :assign_out: Where to save the assignment file
    """
    log.info("Loading trajectories")
    traj_fns = os.listdir(args.traj_dir)
    traj_fns = [os.path.join(args.traj_dir, tfn) for tfn in traj_fns if tfn.endswith('.h5')]
    log.debug("Loading trajectories from %s", str(traj_fns))
    trajs = [mdtraj.load(tfn, stride=args.load_stride) for tfn in traj_fns]

    log.info("Starting cluster")
    metric = toy.Euclidean2d()
    hkm = clustering.KMeans(metric, trajs, distance_cutoff=args.distance_cutoff)

    assignments = clustering.split(hkm._assignments, hkm._traj_lengths)
    assignments = np.array(assignments, dtype=int)
    np.save(args.assign_fn, assignments)


def implied_timescales_func(args):
    """Make implied timescales at a particular lag time

    Properties of args:
        :lagtime: The lag time
        :assign_fn: Assignments filename
        :n_it: number of implied timescales
    """
    assignments = np.load(args.assign_fn)
    counts = msmlib.get_count_matrix_from_assignments(assignments, lag_time=args.lagtime)
    rev_counts, t_matrix, populations, mapping = msmlib.build_msm(counts,
                                                                  ergodic_trimming=True,
                                                                  symmetrize='mle')

    n_eigenvectors = args.n_it + 1
    e_values = msma.get_reversible_eigenvectors(t_matrix, n_eigenvectors, populations=populations)[0]

    # Correct for possible change in n_eigenvectors from trimming
    n_eigenvectors = len(e_values)
    n_implied_times = n_eigenvectors - 1

    # make sure to leave off equilibrium distribution
    lag_times = args.lagtime * np.ones((n_implied_times))
    imp_times = -lag_times / np.log(e_values[1: n_eigenvectors])

    lagimp = np.vstack((lag_times, imp_times)).transpose()
    np.save('lagtime-%d.npy' % args.lagtime, lagimp)


def make_jobs_func(args):
    """Make a bunch of jobs for qsub

    Properties of args:
        :begin:
        :end:
        :increment:
        :make_cluster:
    """
    lagtimes = range(args.begin, args.end, args.increment)
    assert np.all(np.asarray(lagtimes) > 0), 'Positive lagtimes only please'

    # Write clustering job
    if args.make_cluster:
        with open('cluster.job', 'w') as f:
            f.write(CLUSTER_JOB_TEMPLATE.format(**vars(args)))

    # Write lagtimes job
    qsubs = []
    for lt in lagtimes:
        with open('lagtime-%d.job' % lt, 'w') as f:
            f.write(IT_JOB_TEMPLATE.format(lagtime=lt, **vars(args)))
        qsubs.append(QSUB_TEMPLATE.format(lagtime=lt))

    # Write submit script
    with open('submit.sh', 'w') as f:
        f.write(SUBMIT_TEMPLATE.format(qsubs=''.join(qsubs)))

    # Make it executable
    os.chmod('submit.sh', os.stat('submit.sh').st_mode | stat.S_IEXEC)


def plot_func(args):
    """Gather and plot implied timescales.

    Properties of args:
        :xmin:
        :xmax:
    """
    lagtimes = np.zeros((0, 2))

    # Join

    fns = os.listdir('.')
    for fn in fns:
        if re.match('^lagtime-[0-9]+.npy$', fn):
            lagtime = np.load(fn)
            lagtimes = np.append(lagtimes, lagtime, axis=0)

    np.save('lagtime-all.npy', lagtimes)

    # Plot
    pp.scatter(lagtimes[:, 0], lagtimes[:, 1])
    pp.xlabel('lagtime')
    pp.ylabel('implied timescale')

    xmin, xmax = pp.xlim()
    if args.xmin is not None:
        xmin = args.xmin
    if args.xmax is not None:
        xmax = args.xmax
    pp.xlim((xmin, xmax))
    pp.yscale('log')

    fig_fn = 'it-plot-%d-%d.png' % (int(xmin), int(xmax))
    pp.savefig(fig_fn)

def parse():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Get implied timescales',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sp = parser.add_subparsers()

    # -------------------------------------------------------------------------
    # Clustering
    # -------------------------------------------------------------------------
    cluster_p = sp.add_parser('cluster', help='Do clustering',
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cluster_p.add_argument('--load_stride', help='stride trajs by this', default=1, type=int)
    cluster_p.add_argument('--traj_dir', help='directory of trajectories', default='../../trajs/')
    cluster_p.add_argument('--distance_cutoff', help='Distance cutoff for number of clusters',
                           default=0.2, type=float)
    cluster_p.add_argument('--assign_fn', help='Save the assignments here',
                           default='assignments.npy')
    cluster_p.set_defaults(func=cluster_func)

    # -------------------------------------------------------------------------
    # Implied Timescales
    # -------------------------------------------------------------------------
    it_p = sp.add_parser('it', help='Get implied timescales',
                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    it_p.add_argument('--assign_fn', help='Assignments file', default='assignments.npy')
    it_p.add_argument('lagtime', help='The lag time to calculate', type=int)
    it_p.add_argument('--n_it', help='Number of implied timscales', type=int,
                      default=6)
    it_p.set_defaults(func=implied_timescales_func)

    # -------------------------------------------------------------------------
    # Make jobs
    # -------------------------------------------------------------------------
    make_jobs_p = sp.add_parser('make_jobs', help='Make a bunch of jobs for a range of lagtimes',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    make_jobs_p.add_argument('begin', help='first lagtime', type=int)
    make_jobs_p.add_argument('end', help='last lagtime', type=int)
    make_jobs_p.add_argument('-i', '--increment', help='increment', type=int, default=1)
    make_jobs_p.add_argument('-s', '--src_dir', help='source directory', required=True)
    make_jobs_p.add_argument('-c', '--make_cluster', action='store_true',
                             help='Make a cluster job as well if this argument is present')
    make_jobs_p.set_defaults(func=make_jobs_func, make_cluster=False)

    # -------------------------------------------------------------------------
    # Gather and plot
    # -------------------------------------------------------------------------
    plot_p = sp.add_parser('plot', help='Plot implied timescales',
                           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    plot_p.add_argument('-xi', '--xmin', help='Minimum. If not specified, lowest', type=float)
    plot_p.add_argument('-xa', '--xmax', help='Maximum for plot.', type=float)
    plot_p.set_defaults(func=plot_func)



    args = parser.parse_args()
    log.debug(str(args))

    args.func(args)

if __name__ == "__main__":
    log.basicConfig(level=log.INFO)
    parse()
