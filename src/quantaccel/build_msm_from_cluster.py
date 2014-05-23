#!/usr/bin/env python
'''
Created on Nov 18, 2013

@author: harrigan
'''
from __future__ import division

import sys

from msmbuilder import MSMLib as msml, clustering
from quantaccel import toy
from maccelerator import model
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
# - mk7 (for muller3) uses a distance cutoff of 0.1
#==============================================================================

NPOINTS = 50
PERCENTS = np.linspace(0.05, 1.0, NPOINTS)


def do(round_i, how):
    """Build an msm by first clustering.

    :param round_i: The round (used for loading trajectories)
    :param how: Whether to do by percent or round.
    """
    if how == 'rnew':
        wall_steps, trajs = model.load_trajectories(round_i,
                                                    load_func=model.load_muller)
    elif how == 'pnew':
        wall_steps, trajs = model.load_trajectories_percent(PERCENTS[round_i],
                                                            load_func=model.load_muller)

    metric = toy.Euclidean2d()
    # Why is lagtime = 19?
    # This is probably from when I was doing lower lagtimes.
    lag_time = 19
    distance_cutoff = 0.1

    log.info("Starting cluster")
    hkm = clustering.KMeans(metric, trajs, distance_cutoff=distance_cutoff)
    assignments = clustering.split(hkm._assignments, hkm._traj_lengths)
    assignments = np.array(assignments)

    # Get counts
    counts = msml.get_count_matrix_from_assignments(assignments,
                                                    lag_time=lag_time)

    # Make the transition matrix
    try:
        _, t_matrix, _, mapping = msml.build_msm(counts, ergodic_trimming=True,
                                                 symmetrize='mle')
    except:
        # If mle doesn't work, fall back on transpose
        _, t_matrix, _, mapping = msml.build_msm(counts, ergodic_trimming=True,
                                                 symmetrize='transpose')

    # Make a folder for msm's
    try:
        os.mkdir('msms')
    except OSError as e:
        log.debug('Folder exists %s', str(e))

    # Filenames
    centroids_fn = os.path.join('msms',
                                'centroids-%s-mk7-%d.npy' % (how, round_i))
    tmat_fn = os.path.join('msms',
                           'tmatfromclus-%s-mk7-%d.mtx' % (how, round_i))

    log.info("Saving trimmed centroids.")
    trimmed_centroids = hkm._centroids[np.where(mapping != -1)[0]]
    np.savetxt(centroids_fn, trimmed_centroids)

    with open(tmat_fn, 'w') as f:
        scipy.io.mmwrite(f, t_matrix, comment='Wallsteps: %d' % wall_steps)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print """Usage: xxx.py round_i how.
                    how = {round, percent, rnew, pnew}"""

    round_i = int(sys.argv[1])
    how = sys.argv[2]
    do(round_i, how)
