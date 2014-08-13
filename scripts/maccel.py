__author__ = 'harrigan'

import argparse
import maccelerator as maccel
import os
from os.path import join as pjoin
from pkg_resources import resource_filename
import shutil
import glob

import logging

logging.basicConfig(level=logging.INFO)

SHM = '/dev/shm/'


def get_fn(fn):
    return resource_filename('maccelerator', 'reference/{}'.format(fn))


def alanine_entry(args):
    """Entry point for argparse.

    :param args: argparse Namespace or similar
    """
    run_alanine(args.run_id, args.parallel)


def run_alanine(run_id, parallel):
    """Run a grid of alanine configurations

    :param run_id: A label for this copy of a run. Will be used as folder name
    :param parallel: How parallel to do it.
    """
    griddir = 'copy-{}'.format(run_id)
    shm_griddir = pjoin(SHM, griddir)

    # Do it in shared memory
    assert not os.path.exists(griddir)
    os.mkdir(shm_griddir)

    # Set up configuration
    config = maccel.AlanineConfiguration(get_fn('ala.msm.pickl'),
                                         get_fn('ala.centers.h5'))

    grid = maccel.MAccelGrid(config, shm_griddir, run_id, parallel=parallel)
    grid.grid()

    for traj_dir in glob.iglob(pjoin(shm_griddir, '*/trajs/')):
        logging.info('Archiving %s', traj_dir)
        shutil.make_archive(base_name=pjoin(traj_dir, '..', 'trajs'),
                            format='gztar',
                            root_dir=pjoin(traj_dir, '..'), base_dir='trajs')
        shutil.rmtree(traj_dir)

    shutil.copytree(shm_griddir, griddir)
    shutil.rmtree(shm_griddir)


def parse():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--parallel', action='store_true', default=False)

    sp = parser.add_subparsers()

    alanine = sp.add_parser('alanine')
    alanine.add_argument('--run_id', '-i', help='Run ID',
                         type=int, default=0)
    alanine.set_defaults(func=alanine_entry)

    args = parser.parse_args()
    print(args)

    try:
        args.func(args)
    except AttributeError:
        parser.parse_args(['-h'])


if __name__ == "__main__":
    parse()