#!/usr/bin/env python
'''
Created on Nov 18, 2013

@author: harrigan
'''
from __future__ import division

import sys

from msmbuilder import MSMLib as msml, clustering
from quantaccel import toy
from maccelerator import maccel
import scipy.io

import logging as log
import numpy as np
import os


 #==============================================================================
 # - mk4 used a distance cutoff of 0.2
 # - mk5 uses a distance cutoff of 0.4. This will reduce the number of
 #   states, which will increase the bias but potentially reduce statistical
 #   error
 # - mk6 (for muller3) uses a distance cutoff of 0.3 and a lag time of 19
 #==============================================================================

NPOINTS = 50
PERCENTS = np.linspace(0.05, 1.0, NPOINTS)

def do(round_i, how):

    if how == 'round':
        rounds = toy._get_trajs_fns(".")
        wall_steps, trajs = toy.load_trajs_by_round(rounds, round_i)
    elif how == 'percent':
        wall_steps, trajs = toy.load_trajs_by_percent(".", PERCENTS[round_i])
    elif how == 'rnew':
        wall_steps, trajs = maccel.load_trajectories(round_i)

    
    metric = toy.Euclidean2d()
    lag_time = 19
    distance_cutoff = 0.3
    

    log.info("Starting cluster")
    hkm = clustering.KMeans(metric, trajs, distance_cutoff=distance_cutoff)
    assignments = clustering.split(hkm._assignments, hkm._traj_lengths)
    assignments = np.array(assignments)


    counts = msml.get_count_matrix_from_assignments(assignments, lag_time=lag_time)
    
    try:
        _, t_matrix, _, mapping = msml.build_msm(counts, ergodic_trimming=True,
                                                 symmetrize='mle')
    except:
        _, t_matrix, _, mapping = msml.build_msm(counts, ergodic_trimming=True,
                                                 symmetrize='transpose')    


    try:
        os.mkdir('msms')
    except OSError as e:
        log.debug('Folder exists %s', str(e))

    log.info("Saving trimmed centroids.")
    trimmed_centroids = hkm._centroids[np.where(mapping != -1)[0]]
    np.savetxt(os.path.join('msms', 'centroids-%s-mk6-%d.npy' % (how, round_i)), trimmed_centroids)


    with open(os.path.join('msms', 'tmatfromclus-%s-mk6-%d.mtx' % (how, round_i)), 'w') as f:
        scipy.io.mmwrite(f, t_matrix, comment='Wallsteps: %d' % wall_steps)



if __name__ == "__main__":
    if len(sys.argv) != 4:
        print """Usage: xxx.py round_i how.
                    how = {round, percent, rnew}"""

    round_i = int(sys.argv[1])
    how = sys.argv[2]

    do(round_i, how)
