__author__ = 'harrigan'

import unittest
import argparse
from maccelerator.test_utils import get_fn


def parse():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', '-v', action='count', default=1)
    args = parser.parse_args()
    main(args.verbose)


def main(verbosity=1):
    """Discover and run tests."""
    runner = unittest.TextTestRunner(verbosity=verbosity)
    suite = unittest.defaultTestLoader.discover(
        get_fn('../../testing/'))
    runner.run(suite)


if __name__ == "__main__":
    parse()
