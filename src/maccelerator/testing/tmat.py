"""Test things related to transition matrix sampling."""


import unittest
import maccelerator as maccel

class TestTMatSampling(unittest.TestCase):

    def setUp(self):
        self.tmat = None

    def test_sample(self):

        seq = self.tmat.sample(5)

        # Check length
        self.assertEqual(len(seq), 5)

        # TODO: Check sequence


if __name__ == "__main__":
    unittest.main()