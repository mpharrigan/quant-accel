"""Do adaptive sampling with a provided transition matrix (no molecular
dynamics). We can fit the convergence over time to compute a 'speed up' vs.
the parameters in adaptive sampling."""

import scipy.io
import scipy.linalg
import scipy.sparse.linalg
import numpy as np
import logging as log
import pickle
import itertools

from msmbuilder import MSMLib as msmlib
from matplotlib import pyplot as pp


class TMatSimulator(object):

    """Hold the transition matrix and run dynamics."""

    @property
    def n_states(self):
        """Number of states in the model."""
        return self.t_matrix.shape[0]

    def __init__(self, tmat_fn='ntl9.mtx'):
        # Load transition matrix
        t_matrix = scipy.io.mmread(tmat_fn)
        t_matrix = t_matrix.tocsr()
        self.t_matrix = t_matrix

        self.vals, self.vecs = scipy.sparse.linalg.eigs(t_matrix)

        log.info('Loaded transition matrix of shape %s',
                 self.t_matrix.shape)

        # Load generators
        #self.gens = mdtraj.load(self.gens_fn)

    def simulate(self, state_i, number_of_steps=10000, report_interval=1000,
                 out_fn=None):
        """We run some KMC dynamics, and then send back the results.

            state_i - initial state

        """
        log.debug('Starting TMat simulation...')

        t_matrix = self.t_matrix

        state_out = np.zeros(number_of_steps // report_interval, dtype=int)

        report = report_interval
        report_toti = 0

        for _ in xrange(number_of_steps):
            # Get stuff from our sparse matrix

            csr_slicer = slice(
                t_matrix.indptr[state_i],
                t_matrix.indptr[state_i + 1])
            probs = t_matrix.data[csr_slicer]
            colinds = t_matrix.indices[csr_slicer]

            # Check for normalization
            np.testing.assert_almost_equal(np.sum(probs), 1,
                                           err_msg="TMatrix isn't row normalized.")

            # Find our new state and translate to actual indices
            prob_i = np.sum(np.cumsum(probs) < np.random.rand())
            state_i = colinds[prob_i]

            # Check to see if we report
            if report == report_interval:
                state_out[report_toti] = state_i

                # Reset
                report = 0
                report_toti += 1

            report += 1

        # Write
        assert report_toti == state_out.shape[0], "I did my math wrong."

        if out_fn is not None:
            np.savetxt(out_fn, state_out)
        log.debug('Finished TMat simulation.')

        return state_out


class MSM(object):

    """Hold all the trajectories and build and sample from msm."""

    def __init__(self, lag_time=1, beta=0):
        self.tmat = None
        self.counts = None

        self.traj_list = list()

        self.lag_time = lag_time
        self.beta = beta

    def add(self, traj):
        """Add a new trajectory from which the next MSM will be built."""
        self.traj_list.append(traj)

    def build(self, sim):
        """Build an MSM

        This method uses all the trajectories added to the object so far
        and uses `sim` to get the number of states for easier convergence
        calculation.
        """
        counts = msmlib.get_count_matrix_from_assignments(
            np.array(self.traj_list),
            n_states=sim.n_states,
            lag_time=self.lag_time,
            sliding_window=True)

        _, tmat, _, _ = msmlib.build_msm(counts,
                                         symmetrize='Transpose',
                                         ergodic_trimming=False)
        self.counts = counts
        self.tmat = tmat

    def adapt(self, n_new):
        """Return new starting states

            n_new - number of new states
            beta - 'Temperature': \beta = 1 --> even
                                  \beta > 1 --> refinement
                                  \beta < 1 --> exploration
        """

        counts_per_state = np.array(
            self.counts.sum(axis=1)).flatten() + 10. ** -8
        weights = np.power(counts_per_state, self.beta - 1.0)
        weights /= np.sum(weights)

        cum_weights = np.cumsum(weights)
        new_states = -1 * np.ones(n_new, dtype=int)

        for i in range(n_new):
            new_states[i] = np.sum(cum_weights < np.random.rand())

        log.debug("Starting from states: %s", new_states)

        return new_states

    def error_fro(self, sim):
        """Frobenius norm between transition matrices."""
        diff = self.tmat - sim.t_matrix
        return scipy.linalg.norm(diff.todense(), ord='fro')

    def error_kl(self, sim):
        """KL-divergence between 0-th eigenvector (equilibrium distribution).
        """
        vals, vecs = scipy.sparse.linalg.eigs(self.tmat)
        q = vecs[:, 0]
        p = sim.vecs[:, 0]

        q /= np.sum(q)
        p /= np.sum(p)

        return np.sum(np.where(np.abs(p) > 1.e-8, p * np.log(p / q), 0))

    def error(self, sim, func=error_kl):
        """Helper function that calls a specific error (convergence) function.
        """
        return func(self, sim)


class Accelerator(object):

    """Runs rounds of adaptive sampling and keeps the convergence."""

    def __init__(self, simulator, msm, n_rounds=20):
        self.sim = simulator
        self.msm = msm

        self.n_rounds = n_rounds
        self.errors = list()

    def accelerator_loop(self, n_tpr, n_spt):
        """Run the accelerator loop
            n_rounds - number of rounds
            n_tpr - # trajectories per round (paralellization)
            n_spt - # of steps (length) per trajectory
        """
        n_rounds = self.n_rounds
        starting_states = np.random.randint(low=0, high=self.sim.n_states,
                                            size=n_tpr)

        for round_i in range(n_rounds):
            for traj_i in range(n_tpr):
                traj = self.sim.simulate(number_of_steps=n_spt,
                                         state_i=starting_states[traj_i],
                                         report_interval=1)
                self.msm.add(traj)

            walltime = n_spt * (round_i + 1)
            self.msm.build(self.sim)
            self.errors.append([walltime, self.msm.error(self.sim)])
            starting_states = self.msm.adapt(n_tpr)
            log.info("Built model %4d / %4d", (round_i + 1), n_rounds)

        self.errors = np.array(self.errors)


class RunResult(object):

    """Hold the results of an adaptive run for pickling."""

    def __init__(self, params, errors):
        self.params = params
        self.errors = errors

    def plot(self):
        """Plot the convergence vs. 'walltime'."""
        pp.plot(self.errors[:, 0], self.errors[:, 1],
                'o-', label=self.params['n_spt'])


def main(run_i = -1):
    """Define our parameter sets and run the simulations.

    run_i is for pbsdsh. If it is less than 0, all will be run
    """
    tmat_sim = TMatSimulator('ntl9.mtx')

    defaults = {'lag_time': 1, 'beta': 0, 'n_rounds': 20,
                'n_tpr': 20, 'n_spt': 1000}

    beta = [0, 1, 2]

    params = [
        {'n_spt': 2, 'n_rounds': 200},
        {'n_spt': 10, 'n_rounds': 200},
        {'n_spt': 100, 'n_rounds': 200},
        {'n_spt': 1000, 'n_rounds': 20},
        {'n_spt': 10000, 'n_rounds': 5},
        {'n_tpr': 1},
        {'n_tpr': 10},
        {'n_tpr': 100},
        {'n_tpr': 1000}
    ]

    multierrors = list()
    i = 0

    for (setparam, setbeta) in itertools.product(params, beta):

        if run_i < 0 or run_i == i:
            param = dict(defaults)
            param.update(setparam)
            param.update(beta=setbeta)
            msm = MSM(lag_time=param['lag_time'], beta=param['beta'])
            accelerator = Accelerator(tmat_sim, msm, n_rounds=param['n_rounds'])
            accelerator.accelerator_loop(
                n_tpr=param['n_tpr'],
                n_spt=param['n_spt'])

            rr = RunResult(param, accelerator.errors)
            multierrors.append(rr)
            with open('result-a-%d.pickl' % i, 'w') as f:
                pickle.dump(rr, f)

        i += 1



if __name__ == "__main__":
    log.basicConfig(level=log.INFO)
    main()
