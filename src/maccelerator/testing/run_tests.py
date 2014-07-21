__author__ = 'harrigan'

import unittest
from maccelerator.testing.utils import get_fn


if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    suite = unittest.defaultTestLoader.discover(
        get_fn('../src/maccelerator/testing/'))
    runner.run(suite)
