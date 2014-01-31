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

import logging as log


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
        distance_cutoff = 0.2
    elif which == 'tmat':
        metric = rmsd.RMSD()
        raise Exception("No longer supported.")
        lag_time = 1
        distance_cutoff = 0.2


    log.info("Starting cluster")
    hkm = clustering.KMeans(metric, trajs, distance_cutoff=distance_cutoff)
    centroids = hkm._centroids
    assignments = clustering.split(hkm._assignments, hkm._traj_lengths)
    assignments = np.array(assignments)



    counts = msml.get_count_matrix_from_assignments(assignments, lag_time=lag_time)
    _, t_matrix, _, mapping = msml.build_msm(counts, ergodic_trimming=True,
                                             symmetrize='transpose')


    log.info("Saving trimmed centroids.")
    trimmed_centroids = hkm._centroids[np.where(mapping != -1)[0]]
    np.savetxt('centroids-%s-mkiii-%d.npy' % (how, round_i), trimmed_centroids)


    with open('tmatfromclus-%s-mkiii-%d.mtx' % (how, round_i), 'w') as f:
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
