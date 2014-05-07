'''
Created on Nov 18, 2013

@author: harrigan
'''
from __future__ import division
import sys
import numpy as np
import scipy.io
from maccelerator import model
from msmbuilder import MSMLib as msml
import os
import logging as log

NPOINTS = 50
PERCENTS = np.linspace(0.05, 1.0, NPOINTS)


def do(round_i, how):
    """Build an msm from assignments.

    :param round_i: The round (used for loading trajectories)
    :param how: Whether to do by percent or round.
    """

    if how == 'rnew':
        wall_steps, trajs = model.load_trajectories(round_i,
                                                    load_func=model.load_tmat)
    elif how == 'pnew':
        wall_steps, trajs = model.load_trajectories_percent(PERCENTS[round_i],
                                                            load_func=model.load_tmat)

    lag_time = 1

    for traj in trajs:
        if len(traj) <= 0:
            return

    # Get counts
    counts = msml.get_count_matrix_from_assignments(np.array(trajs),
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
    mapping_fn = os.path.join('msms',
                              'mapping-%s-mk7-%d.npy' % (how, round_i))
    tmat_fn = os.path.join('msms',
                           'tmatfromclus-%s-mk7-%d.mtx' % (how, round_i))

    with open(tmat_fn, 'w') as f:
        scipy.io.mmwrite(f, t_matrix, comment='Wallsteps: %d' % wall_steps)

    # Save the results of trimming
    mapping = np.where(mapping != -1)[0]
    np.savetxt(mapping_fn, mapping, fmt="%d")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print """Usage: xxx.py round_i how.
                    how = [rnew, pnew]"""

    round_i = int(sys.argv[1])
    how = sys.argv[2]

    do(round_i, how)
