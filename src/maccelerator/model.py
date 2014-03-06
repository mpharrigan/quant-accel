'''
Created on Mar 5, 2014

@author: harrigan
'''

from quantaccel import toy
import logging as log
from msmbuilder import MSMLib as msml, clustering
import numpy as np
import os
import mdtraj as md



def model(trajs, lagtime, distance_cutoff):
    """Get counts
    
    returns
        :centroids: a trajectory
    """

    metric = toy.Euclidean2d()

    log.info("Starting cluster")
    hkm = clustering.KMeans(metric, trajs, distance_cutoff=distance_cutoff)
    assignments = clustering.split(hkm._assignments, hkm._traj_lengths)
    assignments = np.array(assignments)

    counts = msml.get_count_matrix_from_assignments(assignments, lag_time=lagtime)
    #_, t_matrix, _, mapping = msml.build_msm(counts, ergodic_trimming=True, symmetrize='mle')
    
    centroids = hkm._centroids
    
    # Hack out a trajectory
    centroids_t = trajs[0]
    centroids_t.xyz = centroids[:, np.newaxis, :]
    
    return counts, centroids_t
    
