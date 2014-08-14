__author__ = 'harrigan'

import logging

from maccelerator.testing import run_alanine
import maccelerator as maccel
from os.path import join as pjoin

log = logging.getLogger(__name__)


def do():
    # Do a run, but only get the folder name
    # Simulate starting from just data
    rundir = run_alanine.run()

    run = maccel.MAccelRun.load(pjoin(rundir, 'run.pickl'))

    log.info('Making Movie')
    movie = maccel.PlotMaker(run)
    movie.make_plots()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    do()
