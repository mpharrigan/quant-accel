__author__ = 'harrigan'

import unittest

import maccelerator as maccel
import numpy as np


class TestUniformAdapt(unittest.TestCase):
    def test_uniform(self):
        found_states = np.arange(5) * 2 + 1
        resinds = maccel.UniformAdapter._adapt(found_states, 7)

        self.assertEqual(len(resinds), 7)

        for ri in resinds:
            # Only odd states
            self.assertEqual(ri % 2, 1)


if __name__ == "__main__":
    unittest.main()
