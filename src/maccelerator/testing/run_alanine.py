__author__ = 'harrigan'

import logging

import maccelerator as maccel
from maccelerator.testing.utils import get_fn, get_folder


def run():
    # Do it
    configuration = maccel.AlanineConfiguration(get_fn('ala.msm.pickl'),
                                                get_fn('ala.centers.h5'))
    param = maccel.AlanineParams(spt=40, tpr=1)
    run = maccel.MAccelRun(configuration, param, get_folder('ala'),
                           parallel=False)
    run.run()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    run()
