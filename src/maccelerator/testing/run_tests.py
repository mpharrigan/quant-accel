__author__ = 'harrigan'

import unittest

if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    suite = unittest.defaultTestLoader.discover('.')
    runner.run(suite)
