"""OpenMM simulation device. We connect to server, request a starting
structure, and then propagate.
"""
#############################################################################
# Imports
##############################################################################

import os

from IPython.utils.traitlets import CInt
import mdtraj
import scipy.io

import numpy as np

from ..core.device import Device
from ..core.traitlets import FilePath



    # configurables.
    tmat_fn = FilePath('tProb.mtx', config=True, help='''
        Path to the transition matrix from which to sample.''')

    gens_fn = FilePath('Gens.h5', config=True, help='''
        Path to the generators trajectory.''')

    number_of_steps = CInt(10000, config=True, help='''
        Number of steps of dynamics to do''')

    report_interval = CInt(1000, config=True, help='''
        Interval at which to save positions to a disk, in units of steps''')


    t_matrix = None
    gens = None

    def start(self):
        # Load transition matrix
        t_matrix = scipy.io.mmread(self.tmat_fn)
        t_matrix = t_matrix.tocsr()
        self.t_matrix = t_matrix
        self.log.info('Loaded transition matrix of shape %s',
                      self.t_matrix.shape)

        # Load generators
        self.gens = mdtraj.load(self.gens_fn)

        super(TMatSimulator, self).start()


    def simulate(self, header, content):
        """Main method that is "executed" by the receipt of the
        msg_type == 'simulate' message from the server.

        We run some KMC dynamics, and then send back the results.
        """
        self.log.info('Starting TMat simulation...')

        # Get that starting state path
        starting_state_path = content.starting_state.path

	state_out = np.zeros(number_of_steps // report_interval)

        report = self.report_interval
        report_toti = 0
        for _ in xrange(self.number_of_steps):
            # Get stuff from our sparse matrix
            t_matrix = self.t_matrix
            csr_slicer = slice(t_matrix.indptr[state_i], t_matrix.indptr[state_i + 1])
            probs = t_matrix.data[csr_slicer]
            colinds = t_matrix.indices[csr_slicer]

            # Check for normalization
            np.testing.assert_almost_equal(np.sum(probs), 1,
                                  err_msg="TMatrix isn't row normalized.")

            # Find our new state and translate to actual indices
            prob_i = np.sum(np.cumsum(probs) < np.random.rand())
            state_i = colinds[prob_i]

            # Check to see if we report
            if report == self.report_interval:
		state_out[report_toti] = state_i

                # Reset
                report = 0
                report_toti += 1

            report += 1

        # Write
        assert report_toti == state_out.shape[0], "I did my math wrong."

	#TODO filename
        np.savetxt('simtmat.dat', state_out) 
        self.log.info('Finished TMat simulation.')






###############################################################################
# Utilities
###############################################################################

def random_seed():
    """Get a seed for a random number generator, based on the current platform,
    pid, and wall clock time.

    Returns
    -------
    seed : int
        The seed is a 32-bit int
    """
    import platform
    import time
    import hashlib
    plt = ''.join(platform.uname())
    seed = int(hashlib.md5('%s%s%s' % (plt, os.getpid(), time.time())).hexdigest(), 16)

    return seed % np.iinfo(np.int32).max
