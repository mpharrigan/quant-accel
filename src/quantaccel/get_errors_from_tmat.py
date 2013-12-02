from __future__ import division

import os
import re
import sys

from msmbuilder import msm_analysis as msma
import scipy.io

import numpy as np


NEIGS = 5
EPS = 1.0e-10

def do(which, how, whence):
    """
    which - muller, tmat
    how - round, percent
    whence - tmatfromass, (tmatfromclus?)
    """

    if which == 'muller':
        lag_time = 20
    elif which == 'tmat':
        lag_time = 1

    files = os.listdir('.')
    matchy = '%s-%s-[0-9]+.mtx' % (whence, how)
    i = 0
    its = np.zeros((len(files), NEIGS + 1))
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
            for j, v in enumerate(vals[1:]):
                if np.abs(v - 1.0) < EPS:
                    its[i, j + 1] = 1.0 / EPS
                else:
                    its[i, j + 1] = -lag_time / np.log(v)
            i += 1

    np.savetxt('its-%s.dat' % whence, its)



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
