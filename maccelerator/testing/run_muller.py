__author__ = 'harrigan'

import logging

import maccelerator as maccel
from maccelerator.testing.test_utils import get_fn, get_folder


def run():
    configuration = maccel.MullerConfiguration(get_fn('muller_sys.xml'),
                                               get_fn('muller_int.xml'))
    param = maccel.MullerParams(spt=800, tpr=2)
    run = maccel.MAccelRun(configuration, param, get_folder('mul'))
    run.run()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    run()
