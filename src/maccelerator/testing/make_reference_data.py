__author__ = 'harrigan'

import os
from os.path import join as pjoin
import argparse

from maccelerator.configurations.muller import generate_muller_sysint
from maccelerator.simulate import serialize_openmm


def make_muller_reference_data(dirname):
    system, integrator = generate_muller_sysint()
    serialize_openmm(system, integrator, pjoin(dirname, 'muller_sys.xml'),
                     pjoin(dirname, 'muller_int.xml'))


def make_reference_data(dirname='../../reference'):
    try:
        os.mkdir(dirname)
    except OSError:
        pass

    make_muller_reference_data(dirname)


def parse():
    parser = argparse.ArgumentParser(
        description="Generate reference data for quant accel",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('dirname', help='Directory to write data',
                        default='./reference')

    args = parser.parse_args()
    make_reference_data(args.dirname)


if __name__ == "__main__":
    parse()
