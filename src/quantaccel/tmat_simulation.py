"""Do adaptive sampling with a provided transition matrix (no molecular
dynamics). We can fit the convergence over time to compute a 'speed up' vs.
the parameters in adaptive sampling."""

from __future__ import division

import scipy.io
import scipy.linalg
import scipy.sparse.linalg
import numpy as np
import logging as log
import pickle
import itertools
import argparse
import scipy.optimize

import scipy.sparse.csgraph
import scipy.sparse

from msmbuilder import MSMLib as msmlib
from mdtraj import io

from scipy.sparse.linalg.eigen.arpack import ArpackNoConvergence, ArpackError

STARTING_STATE_FN = 'starting_state.int'
POPCONVERGE = 0.4
NPROCS = 16
NADAPTIVE = 90


class TMatSimulator(object):
    """Hold the transition matrix and run dynamics.


    :param get_eigen: Whether to precompute eigensystem. Defaults to true
                to maintain compatibility with previous system.
    """

    @property
    def n_states(self):
        """Number of states in the model."""
        return self.t_matrix.shape[0]

    def __init__(self, tmat_fn='../ntl9.mtx', get_eigen=True):
        # Load transition matrix
        t_matrix = scipy.io.mmread(tmat_fn)
        t_matrix = t_matrix.tocsr()
        self.t_matrix = t_matrix

        if get_eigen:
            self.p, actual_lam = _get_eigenvec(t_matrix, eigenval=True)
            self.actual_it = -1.0 / np.log(actual_lam)
            self.p /= np.sum(self.p)

        log.info('Loaded transition matrix of shape %s',
                 self.t_matrix.shape)

        # Load generators
        # self.gens = mdtraj.load(self.gens_fn)

    def simulate(self, state_i, number_of_steps, out_fn=None):
        """We run some KMC dynamics, and then send back the results.

        :param state_i: Initial state index
        :param number_of_steps: Length of trajectory to return. Note:
                        We actually take one fewer *step* because we include
                        the initial state in our trajectory
        :param out_fn: Optionally write out a trajectory in HDF5 format
        """
        log.debug('Starting TMat simulation...')

        t_matrix = self.t_matrix

        state_out = np.zeros(number_of_steps, dtype=int)

        # Set first state to initial state!
        state_out[0] = state_i

        for i in xrange(1, number_of_steps):
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

            state_out[i] = state_i

        # Write
        if out_fn is not None:
            io.saveh(out_fn, state_out)
        log.debug('Finished TMat simulation.')
        return state_out


def _super_debug_get_eq_distr(vals, vecs):
    for i, val in enumerate(vals):
        if i > 0:
            log.warn("The first eigenvec is no good. Trying others.")

        if np.abs(val - 1.0) < 1e-8:
            q = vecs[:, i]
            num_pos = len(np.where(q > 1e-8)[0])
            num_neg = len(np.where(q < -1e-8)[0])

            log.debug('Checking vector %d, val %f. Num +/- %d, %d',
                      i, val, num_pos, num_neg)

            if not (num_pos != 0 and num_neg != 0):
                # This is a good one
                break
        else:
            raise ValueError("""We've gotten to an eigenvec with
                                eigenvalue != 1.""")
    else:
        # This is if we make it through the for loop
        # (python has nifty constructs like this)
        raise ValueError("""We've exhausted all eigenvecs.""")

    return q


def _get_eigenvec(tmat, eigenval=False):
    try:
        ttr = tmat.transpose()
        vals, vecs = scipy.sparse.linalg.eigs(ttr, which="LR",
                                              maxiter=10000, tol=1e-30)

        order = np.argsort(-np.real(vals))
        vals = np.real_if_close(vals[order])
        vecs = np.real_if_close(vecs[:, order])

        log.debug('Eigenvalues: %s', str(vals))

        num_unity_eigenvals = np.sum(np.abs(vals - 1.0) < 1e-8)
        log.debug("Number of unity eigenvals: %d", num_unity_eigenvals)
        if num_unity_eigenvals != 1:
            raise ValueError("We found %d eigenvalues that are 1",
                             num_unity_eigenvals)

        q = vecs[:, 0]
    except (ArpackNoConvergence, ArpackError, ValueError) as e:
        log.warn("No eigenv convergence: %s", str(e))
        log.warn("Returning uniform distribution.")
        q = np.ones(ttr.shape[0])
        vals = [1, 1e-8]

    if eigenval:
        return q, vals[1]
    else:
        return q


def _invert_mapping(mapping):
    backmap = -2 * np.ones(np.max(mapping) + 1)
    for bmap, fmap in enumerate(mapping):
        if fmap > -1:
            backmap[fmap] = bmap
    return backmap


class MSM(object):
    """Hold all the trajectories and build and sample from msm."""

    def __init__(self, lag_time=1, beta=0):
        self.tmat = None
        self.counts = None

        self.traj_list = list()

        self.lag_time = lag_time
        self.beta = beta
        self.n_states = 0

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

        try:
            _, tmat, _, mapping = msmlib.build_msm(counts,
                                                   symmetrize='MLE',
                                                   ergodic_trimming=True)
        except Exception as e:
            log.warn("Building msm with MLE failed: %s", str(e))
            log.warn("Defaulting to transpose")
            _, tmat, _, mapping = msmlib.build_msm(counts,
                                                   symmetrize='Transpose',
                                                   ergodic_trimming=False)

        self.counts = counts
        self.n_states = tmat.shape[0]

        # Back out full transition matrix
        backmap = _invert_mapping(mapping)
        coo_tmat = tmat.tocoo()
        back_row = backmap[coo_tmat.row]
        back_col = backmap[coo_tmat.col]
        expanded_tmat = scipy.sparse.coo_matrix(
            (coo_tmat.data, (back_row, back_col)),
            shape=(sim.n_states, sim.n_states))

        self.tmat = expanded_tmat

    def adapt(self, n_new, adapt_func):
        if adapt_func in [MSM.adapt_weights, MSM.adapt_sort]:
            # Do adaptive sampling.
            # Note: my old version of the adaptive code is better
            #       than the one based on ergodic trimming et al
            #       because it trims too much.
            return self.adapt_old(n_new, adapt_func)
        else:
            # Do not do adaptive sampling. Just continue previous trajectories.
            return self.dont_adapt(n_new)

    def dont_adapt(self, n_new):
        """Return the ends of the previous n_new trajectories.

        This effectively continues the previous trajectories.
        """
        # For the last n_new trajectories, pluck the last value
        new_states = [traj[-1] for traj in self.traj_list[-n_new:]]
        return np.array(new_states)

    def adapt_new(self, n_new, adapt_func):
        """Use ergodic trimming and backtrack to sample new.

        This trims too much: specifically states that only have a count 'to'
        and states that only have a count 'from' will not be sampled.
        """
        erg_counts, mapping = msmlib.ergodic_trim(self.counts)

        log.debug("shapes of trimmed vs original: %s %s",
                  str(erg_counts.shape), str(self.counts.shape))

        counts_per_state = np.array(erg_counts.sum(axis=1).flatten() + 1e-8)
        weights = np.power(counts_per_state, self.beta - 1.0)
        weights /= np.sum(weights)

        cum_weights = np.cumsum(weights)
        new_states = -1 * np.ones(n_new, dtype=int)

        for i in range(n_new):
            new_state_erg = np.sum(cum_weights < np.random.rand())
            backtrack = np.where(mapping == new_state_erg)[0]
            assert len(backtrack) == 1, 'Messed up mapping.'
            new_states[i] = backtrack[0]

        log.debug("Starting from states: %s", new_states)

        return new_states

    def adapt_old(self, n_new, adapt_func):
        """Return new starting states

            n_new - number of new states
            beta - 'Temperature': \beta = 1 --> even
                                  \beta > 1 --> refinement
                                  \beta < 1 --> exploration
        """

        num_fro = np.array(self.counts.sum(axis=1), dtype=int).flatten()
        num_to = np.array(self.counts.sum(axis=0), dtype=int).flatten()

        no_fro = np.where(num_fro == 0)[0]
        no_to = np.where(num_to == 0)[0]

        abs_no = np.intersect1d(no_fro, no_to)
        only_no_fro = np.setdiff1d(no_fro, no_to)
        only_no_to = np.setdiff1d(no_to, no_fro)

        log.debug("Undiscovered states: %d", len(abs_no))
        log.debug(
            "Only no fro: %5d\tOnly no to: %5d",
            len(only_no_fro),
            len(only_no_to))

        new_states = adapt_func(self, n_new, abs_no)
        log.debug("Starting from states: %s", new_states)

        return new_states

    def adapt_weights(self, n_new, abs_no):
        """Adapt based on weights and beta.

        :n_new: number of new states to generate
        :abs_no: states we should not sample
        """
        counts_per_state = np.array(
            self.counts.sum(axis=1)).flatten() + 1e-8

        weights = np.power(counts_per_state, self.beta - 1.0)

        # Set prob of choosing 'undiscovered' state to 0
        weights[abs_no] = 0.0

        weights /= np.sum(weights)

        cum_weights = np.cumsum(weights)
        new_states = -1 * np.ones(n_new, dtype=int)

        for i in range(n_new):
            new_states[i] = np.sum(cum_weights < np.random.rand())
            if np.in1d(new_states[i], abs_no):
                log.error("Picking an unfound state")
        log.debug("cps %s", str(counts_per_state[new_states]))
        return new_states


    def adapt_sort(self, n_new, abs_no):
        """Adapt by sorting and selecting in order.

        :n_new: number of new states to generate
        :abs_no: states we should not sample
        """
        counts_per_state = np.array(self.counts.sum(axis=1)).flatten()
        states_to_sample = np.argsort(counts_per_state)
        states_to_sample = np.setdiff1d(states_to_sample, abs_no)
        if len(states_to_sample) > n_new:
            states_to_sample = states_to_sample[:n_new]
        elif len(states_to_sample) < n_new:
            n_copy = n_new // len(states_to_sample)
            n_add = n_new % len(states_to_sample)
            states_to_sample = np.repeat(states_to_sample, n_copy)
            states_to_sample = np.append(states_to_sample,
                                         states_to_sample[:n_add])
        if len(states_to_sample) != n_new:
            log.error("Incorrect number of new states")
        return states_to_sample

    def error_fro(self, sim):
        """Frobenius norm between transition matrices."""
        diff = self.tmat - sim.t_matrix
        return scipy.linalg.norm(diff.todense(), ord='fro')

    def error_kl(self, sim):
        """KL-divergence between 0-th eigenvector (equilibrium distribution).
        """
        p = sim.p
        q = _get_eigenvec(self.tmat)

        q /= np.sum(q)
        p /= np.sum(p)

        return np.sum(
            np.nan_to_num(
                np.where(np.abs(p) > 1.e-6, p * np.log(p / q), 0)))

    def error_tvd(self, sim):
        """Total variation distance."""
        p = sim.p
        q = _get_eigenvec(self.tmat)

        norm_to = 1000.0
        q = norm_to * (q / np.sum(q))
        p = norm_to * (p / np.sum(p))

        res = 0.5 * (1 / norm_to) * np.sum(np.abs(p - q))

        return res

    def error_it_logdiff(self, sim):
        """Calculate our implied timescale and subtract in log-space from
        actual one."""
        actual_it = sim.actual_it
        _, est_lam = _get_eigenvec(self.tmat, eigenval=True)
        est_it = -1.0 / np.log(est_lam)
        return np.log(actual_it / est_it)

    def error_euc(self, sim):
        """Euclidean distance."""
        p = sim.p
        q = _get_eigenvec(self.tmat)

        q /= np.sum(q)
        p /= np.sum(p)

        diff = p - q

        return np.sqrt(np.dot(diff, diff))

    def error(self, sim, func=error_tvd):
        """Helper function that calls a specific error (convergence) function.
        """
        return func(self, sim)


class Accelerator(object):
    """Runs rounds of adaptive sampling and keeps the convergence."""

    def __init__(self, simulator, msm):
        self.sim = simulator
        self.msm = msm

        self.poperrors = list()
        self.iterrors = list()
        self.nstates = list()

    def accelerator_loop(self, n_tpr, n_spt, adapt_func):
        """Run the accelerator loop
            n_tpr - # trajectories per round (paralellization)
            n_spt - # of steps (length) per trajectory
            adaptive - whether to run adaptive sampling or not
            :adapt_func: function for adaptive sampling
        """

        with open(STARTING_STATE_FN) as f:
            starting_state = int(f.read())

        # Start everything from one starting state
        starting_states = np.ones(n_tpr, dtype='int') * starting_state

        log.debug("Starting states: %s", str(starting_states))

        # We will run until convergence, then start decrementing this
        # value. rounds_left is the number of additional rounds to do
        # after convergence
        rounds_left = 10
        converged = False
        round_i = 0

        while rounds_left > 0:

            # Simulate
            for traj_i in range(n_tpr):
                traj = self.sim.simulate(n_steps=n_spt,
                                         sstate=starting_states[traj_i])
                self.msm.add(traj)

            # Build MSM
            walltime = n_spt * (round_i + 1)
            self.msm.build(self.sim)

            # Record interesting info
            poperror = self.msm.error_tvd(self.sim)
            iterror = self.msm.error_it_logdiff(self.sim)
            nstate = self.msm.n_states
            self.poperrors.append([walltime, poperror])
            self.iterrors.append([walltime, iterror])
            self.nstates.append([walltime, nstate])

            if not converged and poperror < POPCONVERGE:
                log.info("Convergence achieved at round %4d", (round_i + 1))
                converged = True
                rounds_left -= 1
            elif converged:
                rounds_left -= 1

            # Adapt
            starting_states = self.msm.adapt(n_tpr, adapt_func)
            log.info("Built and sampled model %4d", (round_i + 1))

            round_i += 1

        self.poperrors = np.array(self.poperrors)
        self.iterrors = np.array(self.iterrors)
        self.nstates = np.array(self.nstates, dtype=int)


class RunResult(object):
    """Hold the results of an adaptive run for pickling."""

    def __init__(self, params, accelerator=None):
        self.params = params
        if accelerator is not None:
            self.poperrors = accelerator.poperrors
            self.iterrors = accelerator.iterrors
            self.nstates = accelerator.nstates


def simulate(tmat_sim, defaults, set_beta, set_spt, set_tpr, adapt_func):
    """Run one simulation given the params
     """

    log.info("Setting beta = %f spt = %d tpr = %d", set_beta, set_spt, set_tpr)
    log.info("Using %s as adapt method", adapt_func.func_name)

    # Make param dict from defaults and set our set values
    param = dict(defaults)
    param.update(n_spt=set_spt)
    param.update(n_tpr=set_tpr)
    param.update(beta=set_beta)

    # Make MSM container object
    msm = MSM(lag_time=param['lag_time'], beta=param['beta'])

    # Make accelerator object
    accelerator = Accelerator(tmat_sim, msm)

    # Run it
    accelerator.accelerator_loop(n_tpr=param['n_tpr'],
                                 n_spt=param['n_spt'],
                                 adapt_func=adapt_func)

    # Get the results
    rr = RunResult(param, accelerator)
    return rr


def get_params_adaptive(beta):
    spt = [2, 4, 8, 16, 32, 64]
    tpr = [1, 10, 100, 500, 1000]
    log.info("Number of permutations = %d", len(beta) * len(spt) * len(tpr))
    return itertools.product(beta, spt, tpr)


def get_params_adapt_weights():
    beta = [-1, 1, 3]
    return get_params_adaptive(beta)


def get_params_adapt_sort():
    beta = [0]  # Can be anything
    return get_params_adaptive(beta)


def get_params_control():
    beta = 0  # Can be anything

    # (beta, n_spt, n_tpr)
    params = [
        (beta, 100, 1),
        (beta, 10, 10),
        (beta, 10, 100),
        (beta, 4, 500),
        (beta, 4, 1000)
    ]
    return params


def main(run_i=-1, runcopy=0, adapt_func_str='weights'):
    """Define our parameter sets and run the simulations.

    :run_i: run this permutation. set negative to do all
    :runcopy: append this id to distinguish between copies
    :adapt_func_str: string to choose function for adapting
    """
    # Load up the transition matrix
    tmat_sim = TMatSimulator('../ntl9.mtx')

    # Default parameters, i.e. those that do not vary over the different sims
    defaults = {'lag_time': 1, 'runcopy': runcopy}

    # Beta only makes sense for weights
    if adapt_func_str == 'weights':
        params = get_params_adapt_weights()
        adapt_func = MSM.adapt_weights
    elif adapt_func_str == 'sort':
        params = get_params_adapt_sort()
        adapt_func = MSM.adapt_sort
    elif adapt_func_str == 'non':
        params = get_params_control()
        adapt_func = MSM.dont_adapt
    else:
        log.error("Unrecognized adapt func: %s", adapt_func_str)

    for i, (set_beta, set_spt, set_tpr) in enumerate(params):
        if run_i < 0 or run_i == i:
            rr = simulate(tmat_sim, defaults, set_beta, set_spt, set_tpr,
                          adapt_func=adapt_func)
            with open('result-{adapt}-{runcopy}-{i}.pickl'.format(
                    adapt=adapt_func_str, runcopy=runcopy, i=i), 'w') as f:
                pickle.dump(rr, f, protocol=2)


def parse():
    """Parse command line args."""
    parser = argparse.ArgumentParser(
        description='Perform accelerated sampling on transition matrix',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('run_i', type=int,
                        help="Which permute to run")
    parser.add_argument('--runcopy', type=int,
                        help="ID for this copy", default=0)
    parser.add_argument('--adapt', choices=['weights', 'sort', 'non'],
                        help="What method to use for adaptive",
                        default='weights')
    args = parser.parse_args()

    main(args.run_i, runcopy=args.runcopy, adapt_func_str=args.adapt)


if __name__ == "__main__":
    log.basicConfig(level=log.INFO)
    np.seterr(under='warn')
    parse()
