'''
Created on Nov 18, 2013

@author: harrigan
'''
from __future__ import division

import sys

from msmbuilder import MSMLib as msml, clustering
from msmbuilder.metrics import rmsd
import scipy.io

import numpy as np
from quantaccel import toy


NMED_ITERS = 10

NPOINTS = 20
PERCENTS = np.linspace(0.05, 1.0, NPOINTS)

def do(round_i, which, how):

    if how == 'round':
        rounds = toy._get_trajs_fns(".")
        wall_steps, trajs = toy.load_trajs_by_round(rounds, round_i)
    elif how == 'percent':
        wall_steps, trajs = toy.load_trajs_by_percent(".", PERCENTS[round_i])

    if which == 'muller':
        metric = toy.Euclidean2d()
        lag_time = 20
        distance_cutoff = 0.3
    elif which == 'tmat':
        metric = rmsd.RMSD()
        lag_time = 1
        distance_cutoff = 0.2

    shim_trajs = [toy.ShimTrajectory(traj.xyz) for traj in trajs]
    print "Starting cluster"
    hkm = clustering.HybridKMedoids(metric, shim_trajs, k=None,
                                    distance_cutoff=distance_cutoff,
                                    local_num_iters=NMED_ITERS)
    assignments = hkm.get_assignments()

    if round_i == NPOINTS - 1:
        print "Saving generators"
        gens = hkm.get_generators_as_traj()
        gens = metric.prepare_trajectory(gens)
        np.savetxt('gens.npy', gens)

    counts = msml.get_count_matrix_from_assignments(assignments,
                                                    lag_time=lag_time)
    _, t_matrix, _, _ = msml.build_msm(counts, ergodic_trimming=True,
                                       symmetrize='transpose')


    with open('tmatfromclus-%s-%d.mtx' % (how, round_i), 'w') as f:
        scipy.io.mmwrite(f, t_matrix, comment='Wallsteps: %d' % wall_steps)



if __name__ == "__main__":
    if len(sys.argv) != 4:
        print """Usage: xxx.py round_i which how.
                    which = {muller, tmat}
                    how = {round, percent}"""

    round_i = int(sys.argv[1])
    which = sys.argv[2]
    how = sys.argv[3]

    do(round_i, which, how)
