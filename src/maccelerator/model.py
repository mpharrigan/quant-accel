'''
Created on Mar 5, 2014

@author: harrigan
'''

import os

from msmbuilder import MSMLib as msml, clustering
from quantaccel import toy

import logging as log
import mdtraj as md
import numpy as np


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

    counts = msml.get_count_matrix_from_assignments(
        assignments,
        lag_time=lagtime)
    #_, t_matrix, _, mapping = msml.build_msm(counts, ergodic_trimming=True, symmetrize='mle')

    centroids = hkm._centroids
    centroids_t = md.Trajectory(centroids[:, np.newaxis, :], trajs[0].topology)

    return counts, centroids_t
