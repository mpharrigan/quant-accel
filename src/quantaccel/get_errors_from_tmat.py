from __future__ import division

import os
import re
import sys

from msmbuilder import msm_analysis as msma
import scipy.io

import numpy as np

from quantaccel.tmat_simulation import RunResult
import pickle

GVECS_FN = "../../gold/gold/gold_vecs.pickl"

def errors_kl(vecs, gold_vecs):
    """KL-divergence between 0-th eigenvector (equilibrium distribution).
    """

    q = vecs[:, 0]
    p = gold_vecs[:, 0]

    q /= np.sum(q)
    p /= np.sum(p)

    return np.nan_to_num(np.sum(np.where(np.abs(p) > 1.e-8, p * np.log(p / q), 0)))


NEIGS = 4
EPS = 1.0e-10

def do(which, how, whence):
    """
    which - muller, tmat
    how - round, percent
    whence - tmatfromass, tmatfromclus
    """
    with open(GVECS_FN) as f:
        gold_vecs = pickle.load(f)


    if which == 'muller':
        lag_time = 20
    elif which == 'tmat':
        lag_time = 1

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
                print "Number of wallsteps: %d" % wallsteps

            vals, vecs = msma.get_eigenvectors(tmat, n_eigs=NEIGS + 1)
            its[i, 0] = wallsteps
            its[i, 1] = errors_kl(vecs, gold_vecs)
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
