from __future__ import division

import os
import re
import sys

#from msmbuilder import msm_analysis as msma
import scipy.sparse.linalg
import scipy.io

import numpy as np

from quantaccel.tmat_simulation import RunResult
import pickle
from scipy.sparse.linalg.eigen.arpack import ArpackNoConvergence

import logging as log

GVECS_FN = "../../gold/gold/gold_vecs2.pickl"
GVECS_TMAT_FN = "../../gold-vecs2.pickl"

def error_kl(tmat, gold_vecs):
    """KL-divergence between 0-th eigenvector (equilibrium distribution).
    """
    p = gold_vecs[:,0]
    try:
        vals, vecs = scipy.sparse.linalg.eigs(tmat)
        vecs = np.real_if_close(vecs)
        q = vecs[:,0]
    except ArpackNoConvergence:
        log.warn("No eigenv convergence")
        q = np.ones(p.shape)

    q /= np.sum(q)
    p /= np.sum(p)

    return np.sum(np.nan_to_num(np.where(np.abs(p) > 1.e-6, p * np.log(p / q), 0)))

def errors_tvd(tmat, gold_vecs):
    """Total variation distance."""
    p = gold_vecs[:,0]
    try:
        vals, vecs = scipy.sparse.linalg.eigs(tmat)
        vecs = np.real_if_close(vecs)
        q = vecs[:,0]
    except ArpackNoConvergence:
        log.warn("No eigenv convergence")
        q = np.ones(p.shape)

    q /= np.sum(q)
    p /= np.sum(p)

    return 0.5 * np.sum(np.abs(p-q))

def errors(tmat, gold_vecs, func=errors_tvd):
    return func(tmat, gold_vecs)


NEIGS = 4
EPS = 1.0e-10

def do(which, how, whence):
    """
    which - muller, tmat
    how - round, percent
    whence - tmatfromass, tmatfromclus
    """

    if which == 'muller':
        lag_time = 20
        with open(GVECS_FN) as f:
            gold_vecs = pickle.load(f)
            gold_vecs = np.real_if_close(gold_vecs)
    elif which == 'tmat':
        lag_time = 1
        with open(GVECS_TMAT_FN) as f:
            gold_vecs = pickle.load(f)
            gold_vecs = np.real_if_close(gold_vecs)

    files = os.listdir('.')
    matchy = '%s-%s-[0-9]+.mtx' % (whence, how)
    i = 0

    its = np.zeros((len(files), 2))
    for fn in files:
        if re.match(matchy, fn):
            with open(fn) as f:
                tmat = scipy.io.mmread(f)
                f.seek(0)  # go to comment line
                f.readline()
                try:
                    wallsteps = int(f.readline().strip().split()[-1])
                except ValueError:
                    wallsteps = -1
                log.debug("Number of wallsteps: %d", wallsteps)

            its[i, 0] = wallsteps
            its[i, 1] = errors(tmat, gold_vecs)
            i += 1

    # Fix format of errors array
    its = its[:i]
    its = its.transpose()

    # Hacky get params
    cwd = os.path.abspath(os.curdir)
    folder_name = cwd[cwd.rfind('/')+1:]
    fn_splits = folder_name.split('-')
    params = {fn_splits[0]: fn_splits[1]}
    rr = RunResult(params, its)
    with open("%s-%s.kl.pickl" % (whence,how), 'w') as f:
        pickle.dump(rr, f)

    return its



if __name__ == "__main__":
    print """Usage: xxx.py which how whence.
                which = [muller, tmat]
                how = [round, percent]
                whence = [tmatfromass, tmatfromclus]
                """
    which = sys.argv[1]
    how = sys.argv[2]
    whence = sys.argv[3]

    do(which, how, whence)
