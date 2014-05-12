"""Test things related to transition matrix sampling."""

import unittest
import maccelerator as maccel
import numpy as np
import scipy.sparse
import numpy.testing


class TestTMatSampling(unittest.TestCase):
    def setUp(self):
        tmat = scipy.sparse.csr_matrix(
            np.diag(np.ones(5), 1)  # 6 state linear transition model
        )

        self.simulator = maccel.TMatSimulator(tmat)

    def test_sample(self):
        seq = self.simulator.simulate(sstate=0, n_steps=6)

        # Check length
        self.assertEqual(len(seq), 6)

        # Check states
        np.testing.assert_array_equal(seq, range(6))



if __name__ == "__main__":
    unittest.main()