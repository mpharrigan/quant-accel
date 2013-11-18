'''
Created on Nov 18, 2013

@author: harrigan
'''

from quantaccel import toy
from msmbuilder.metrics import rmsd
import sys
import numpy as np
import mdtraj
import scipy.io

MULLER_GENS = '../../cluster_gens.npy'
TMAT_GENS = '../../Gens.h5'

def do(round_i, which):
    rounds = toy._get_trajs_fns(".")
    wall_steps, trajs = toy.load_trajs_by_round(rounds, round_i)

    if which == 'muller':
        metric = toy.Euclidean2d()
        lag_time = 20
        gens = np.load(MULLER_GENS)
        generators = toy.ShimTrajectory(gens)
    elif which == 'tmat':
        metric = rmsd.RMSD()
        lag_time = 1
        gens = mdtraj.load(TMAT_GENS)
        generators = toy.ShimTrajectory(gens.xyz)

    t_matrix = toy.build_msm(trajs, generators, lag_time, metric)
    with open('tmat-%d.mtx' % round_i, 'w') as f:
        scipy.io.mmwrite(f, t_matrix, comment='Wallsteps: %d' % wall_steps)



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print """Usage: xxx.py round_i which.
                    which = {muller, tmat}"""
    round_i = int(sys.argv[1])
    which = sys.argv[2]
    do(round_i, which)
