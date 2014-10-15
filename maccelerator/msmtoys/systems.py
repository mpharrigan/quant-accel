__author__ = 'harrigan'

import numpy as np
import scipy.sparse

from .systems_baseclasses import TransitionSystem
from .muller import MullerForce


EPSILON = 1e-6


class TmatFromPotential(TransitionSystem):
    def __init__(self):
        super().__init__()
        self.grid = None


    def calculate_transition_matrix(self, potential, resolution, beta,
                                    bounds=None):
        """Calculate a transition matrix from a potential."""
        grid = potential.get_grid(resolution, bounds)

        n_states = grid.shape[1] * grid.shape[2]

        # The transition matrix will be constructed as a sparse COO matrix
        t_matrix_rows = np.zeros((0,))
        t_matrix_cols = np.zeros((0,))
        t_matrix_data = np.zeros((0,))

        neighbors = _neighbors()
        potential = potential.potential(grid[0], grid[1])

        for row_i in range(grid.shape[1]):
            for col_i in range(grid.shape[2]):
                # Loop through each starting point
                pot = potential[row_i, col_i]
                from_state = _state_id(row_i, col_i, grid.shape)
                normalization = 0.0
                # Only do nearest-neighbor
                for neighs in neighbors:
                    nrow_i = row_i + neighs[0]
                    ncol_i = col_i + neighs[1]
                    if (nrow_i < 0 or nrow_i >= grid.shape[1]
                        or ncol_i < 0 or ncol_i >= grid.shape[2]):
                        # Transition probability to states outside our
                        # state space is zero.
                        # This is not in the transition matrix
                        pass
                    else:
                        to_state = _state_id(nrow_i, ncol_i, grid.shape)
                        delta_pot = potential[nrow_i, ncol_i] - pot
                        t_prob = np.exp(-beta * delta_pot)

                        # Store info for our sparse matrix
                        t_matrix_rows = np.append(t_matrix_rows, from_state)
                        t_matrix_cols = np.append(t_matrix_cols, to_state)
                        t_matrix_data = np.append(t_matrix_data, t_prob)

                        normalization += t_prob

        t_matrix = scipy.sparse.coo_matrix(
            (t_matrix_data, (t_matrix_rows, t_matrix_cols)),
            shape=(n_states, n_states)).tocsr()

        # Normalize
        for trow_i in range(n_states):
            rfrom = t_matrix.indptr[trow_i]
            rto = t_matrix.indptr[trow_i + 1]
            normalization = np.sum(t_matrix.data[rfrom:rto])
            if normalization < EPSILON:
                print("No transitions from %d" % (trow_i))
            else:
                t_matrix.data[rfrom:rto] = \
                    t_matrix.data[rfrom:rto] / normalization

        self.grid = grid
        self.tmat = t_matrix


class MullerTmat(TmatFromPotential):
    def __init__(self, resolution, beta):
        super(TmatFromPotential, self).__init__()
        self.calculate_transition_matrix(MullerForce, resolution, beta)


def _neighbors():
    """Get an array of x,y offsets to index nearest neighbors."""
    neighbors = np.zeros((0, 2), dtype='int')
    for row_i in range(-1, 2):
        for col_i in range(-1, 2):
            if not (row_i == 0 and col_i == 0):
                neighbors = np.append(neighbors, [[row_i, col_i]], axis=0)

    return neighbors


def _state_id(x, y, shape):
    """Take a row, column description into a state id."""
    return shape[1] * y + x
