'''
Created on Nov 11, 2013

@author: harrigan
'''

import os, re
import numpy as np
import mdtraj
from msmbuilder import MSMLib as msml, msm_analysis as msma
from msmbuilder import clustering
from msmbuilder.metrics.baseclasses import Vectorized
from msmaccelerator.core import markovstatemodel
import sys
import pickle
from mpi4py import MPI
import sqlite3 as sql


COMM = MPI.COMM_WORLD
MPIRANK = COMM.Get_rank()
MPISIZE = COMM.Get_size()

GOLD_DIR = 'quant/gold'
LL_DIR = 'quant/ll'
LPT_DIR = 'quant/lpt'
CLUSTER_FN = 'cluster_gens.npy'
DISTANCE_CUTOFF = 0.25
MEDOID_ITERS = 0
LOAD_STRIDE = 1
LAG_TIME = 20

def cluster():
    traj_fns = os.listdir(os.path.join(GOLD_DIR, 'trajs/'))
    traj_fns = [os.path.join(GOLD_DIR, 'trajs/', tfn) for tfn in traj_fns]
    print "Loading files."
    trajs = mdtraj.load(traj_fns, stride=LOAD_STRIDE)
    shim_trajs = [ShimTrajectory(traj.xyz) for traj in trajs]
    metric = Euclidean2d()

    print "Starting cluster"
    hkm = clustering.HybridKMedoids(metric, shim_trajs, k=None, distance_cutoff=DISTANCE_CUTOFF, local_num_iters=MEDOID_ITERS)
    generators = hkm.get_generators_as_traj()
    n_clusters = generators['XYZList'].shape[0]
    print "Clustered data into {:,} clusters".format(n_clusters)
    genxyz = generators['XYZList']
    np.save(CLUSTER_FN, genxyz)

def build_msm(trajs, generators, lag_time):
    # Do assignment
    shim_trajs = [ShimTrajectory(traj.xyz) for traj in trajs]
    metric = Euclidean2d()

    # Allocate array
    n_trajs = len(trajs)
    max_traj_len = max([traj.xyz.shape[0] for traj in trajs])
    assignments = -1 * np.ones((n_trajs, max_traj_len), dtype='int')

    # Prepare generators
    pgens = metric.prepare_trajectory(generators)

    for i, traj in enumerate(shim_trajs):
        ptraj = metric.prepare_trajectory(traj)

        for j in xrange(len(traj)):
            d = metric.one_to_all(ptraj, pgens, j)
            assignments[i, j] = np.argmin(d)

    counts = msml.get_count_matrix_from_assignments(assignments, n_states=len(generators), lag_time=lag_time)
    _, t_matrix, _, _ = msml.build_msm(counts, ergodic_trimming=False, symmetrize='transpose')
    return t_matrix

def calculate_errors(n_eigen, comp_tmatrices, wall_steps, gold_tmatrix):

    _, gvecs = msma.get_eigenvectors(gold_tmatrix, n_eigs=n_eigen)
    errors = np.zeros((len(comp_tmatrices), 2))
    errors[:, 0] = wall_steps

    for i in xrange(len(comp_tmatrices)):
        _, vecs = msma.get_eigenvectors(comp_tmatrices[i], n_eigs=n_eigen)
        if gvecs.shape[0] != vecs.shape[0]:
            print "Error: Vectors are not the same shape!"

        error = 0.0
        for j in xrange(n_eigen):
            diff = vecs[:, j] - gvecs[:, j]
            error += np.dot(diff, diff)

        errors[i, 1] = error

    return errors

def load_trajs_by_percent(directory, load_up_to_this_percent):
    """Load trajectories by percentage."""

    assert load_up_to_this_percent <= 1.0, 'Must load less than 100%'

    traj_fns = os.listdir(os.path.join(directory, 'trajs/'))
    traj_fns = [os.path.join(directory, 'trajs/', tfn) for tfn in traj_fns]
    trajs = [mdtraj.load(traj_fn) for traj_fn in traj_fns]

    max_end = trajs[0].shape[0]
    load_end = int(load_up_to_this_percent * max_end)
    print "Loading trajectory up to index {:,} out of {:,}".format(load_end, max_end)
    partial_trajs = [traj[:load_end] for traj in trajs]

    wall_steps = load_end
    return (wall_steps, partial_trajs)

def load_trajs_by_round(rounds, round_i):
    """Load trajectories up to round_i."""

    assert round_i < len(rounds), 'Round index out of range %d' % round_i

    print "Using trajectories after round {}".format(round_i + 1)

    traj_fns = rounds[round_i]
    trajs = [mdtraj.load(traj_fn, stride=LOAD_STRIDE) for traj_fn in traj_fns]

    # Stats
    traj_len = trajs[0].xyz.shape[0]
    n_trajs = len(trajs)
    print "Loaded {} trajs x {:,} length = {:,} frames".format(n_trajs, traj_len, n_trajs * traj_len)
    wall_steps = traj_len * (round_i + 1)
    return wall_steps, trajs

def _get_trajs_fns(directory, db_fn='db.sqlite'):
    """Save trajectory filenames per round of simulation.

    This function loads the msmaccelerator database file
    and uses the models table to enumerate rounds. From each
    model, the trajectories generated up to that point
    are saved.
    """
    connection = sql.connect(os.path.join(directory, db_fn))
    cursor = connection.cursor()
    cursor.execute('select path from models order by time asc')
    model_fns = cursor.fetchall()

    rounds = list()
    for fn in model_fns:
        # Unpack tuple
        fn = fn[0]
        # Get relative path
        fn = fn[fn.rfind('/'):]
        fn = "%s/models/%s" % (directory, fn)
        # Load model
        model = markovstatemodel.MarkovStateModel.load(fn)
        # Get relative paths
        traj_filenames = ["%s/trajs/%s" % (directory, tfn[tfn.rfind('/'):]) for tfn in model.traj_filenames]
        # Save trajectory filenames
        rounds.append(traj_filenames)
        # Cleanup
        model.close()

    return rounds

def get_lpt_dirs():
    lpt_dirs = ['lpt-250-120', 'lpt-500-60', 'lpt-1000-30',
                'lpt-2000-15', 'lpt-6000-5', 'lpt-30000-1']

    lpt_dirs = [os.path.join(LPT_DIR, lpt) for lpt in lpt_dirs]
    return lpt_dirs

def get_ll_dirs():
    ll_dirs = ['ll-2', 'll-5', 'll-6', 'll-7', 'll-10', 'll-14',
               'll-15', 'll-30', 'll-35', 'll-42', 'll-70']
    ll_dirs = [os.path.join(LL_DIR, ll) for ll in ll_dirs]
    return ll_dirs

def lpt_dir_to_x(fn):
    fn = fn[fn.rfind('/') + 1:]
    reres = re.search('(?<=-)[0-9]+(?=-)', fn)
    x = int(reres.group(0))
    return x

def ll_dir_to_x(fn):
    fn = fn[fn.rfind('/') + 1:]
    reres = re.search('(?<=-)[0-9]+', fn)
    x = int(reres.group(0))
    return x


class Euclidean2d(Vectorized):
    def prepare_trajectory(self, trajectory):
        xyz = trajectory["XYZList"]
        if xyz.shape[2] == 3:
            xyz[:, :, 2] = 0.0
        nframes, natoms, ndims = xyz.shape
        return xyz.reshape(nframes, natoms * ndims).astype(np.double)

class ShimTrajectory(dict):
    """This is a dict that can be used to interface some xyz coordinates
    with MSMBuilder's clustering algorithms.

    I'm really sorry that this is necessary. It's horridly ugly, but it comes
    from the fact that I want to use the mdtraj trajectory object (its better),
    but the msmbuilder code hasn't been rewritted to use the mdtraj trajectory
    yet. Soon, we will move mdtraj into msmbuilder, and this won't be necessary.
    """
    def __init__(self, xyz):
        super(ShimTrajectory, self).__init__()
        self['XYZList'] = xyz

    def __len__(self):
        return len(self['XYZList'])

def do_cluster():
    print "Clustering"
    cluster()

def actually_do(directory, generators):
    print "Doing", directory
    rounds = _get_trajs_fns(directory)
    tmatrices = list()
    wall_steps = np.zeros(len(rounds))
    for i in xrange(len(rounds)):
        ws, trajs = load_trajs_by_round(rounds, i)
        tmatrices.append(build_msm(trajs, generators, lag_time=LAG_TIME))
        wall_steps[i] = ws
        
    import pdb; pdb.set_trace()

    with open(os.path.join(directory, 'tmats.pickl'), 'w') as f:
        pickle.dump((directory, wall_steps, tmatrices), f, protocol=2)


def do_tmatrices(which):

    print "Rank %d/%d reporting in!" % (MPIRANK, MPISIZE)

    # Divvy
    if which == 1:
        xx_dirs = get_lpt_dirs()
    elif which == 2:
        xx_dirs = get_ll_dirs()
    else:
        return

    n_tmats = len(xx_dirs)
    n_jobs = n_tmats // MPISIZE
    n_leftovers = n_tmats % MPISIZE

    gens = np.load(CLUSTER_FN)
    generators = ShimTrajectory(gens)

    # Do normal
    for i in xrange(n_jobs):
        dirs = xx_dirs[i * MPISIZE + MPIRANK]
        actually_do(dirs, generators)

    if MPIRANK < n_leftovers:
        dirs = xx_dirs[n_jobs * MPISIZE + MPIRANK]
        actually_do(dirs, generators)
        
def load_tmatrices(which):
    
    with open('quant/lpt/lpt-2000-15/tmats.pickl') as f:
        dir, ws, tmats = pickle.load(f)
    
    print "dir", dir
    print "ws", ws
    print "len", len(tmats)   


if __name__ == "__main__":
    if len(sys.argv) == 2:
        if sys.argv[1] == 'cluster':
            do_cluster()
        elif sys.argv[1] == 'tmatrices':
            do_tmatrices(1)
            do_tmatrices(2)
        elif sys.argv[1] == 'test':
            do_tmatrices(1)
        elif sys.argv[1] == 'errors':
            load_tmatrices(1)
        else:
            print "Not specified"
    else:
        print "Invalid arguments"
