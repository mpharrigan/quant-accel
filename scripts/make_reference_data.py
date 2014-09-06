"""A script to generate reference data for testing.

This can be run from the command line.
"""

import os
from os.path import join as pjoin
import argparse
import pickle

import scipy.io
import mdtraj.io
from mixtape.datasets.alanine_dipeptide import fetch_alanine_dipeptide

from maccelerator.configurations.muller import generate_muller_sysint
from maccelerator.configurations.alanine import generate_alanine_msm
from maccelerator.simulate import serialize_openmm


def make_alanine_reference_data(dirname):
    """Make a small transition matrix from Alanine trajectories

    We featurize using phi / psi angles

    Note: This function is not-deterministic, although it would be useful
    if it were, so testing could be conducted.
    """
    ala_trajs_dir = pjoin(dirname, 'ala_trajs')
    try:
        os.mkdir(ala_trajs_dir)
    except OSError:
        pass

    ala = fetch_alanine_dipeptide(ala_trajs_dir)
    msm, kmeans = generate_alanine_msm(ala)

    # Save cluster centers
    mdtraj.io.saveh(pjoin(dirname, 'ala.centers.h5'),
                    cluster_centers=kmeans.cluster_centers_)

    # Save transition matrix
    scipy.io.mmwrite(pjoin(dirname, 'ala.mtx'), msm.transmat_,
                     comment='Generated for quant-accel reference data')

    # Save MSM Object
    with open(pjoin(dirname, 'ala.msm.pickl'), 'wb') as f:
        pickle.dump(msm, f)


def make_muller_reference_data(dirname):
    """Make OpenMM system.xml and integrator.xml files.

    :param dirname: Where to put the files.
    """

    # TODO: Integrator has a random seed
    system, integrator = generate_muller_sysint()
    serialize_openmm(system, integrator, pjoin(dirname, 'muller_sys.xml'),
                     pjoin(dirname, 'muller_int.xml'))


def make_reference_data(dirname='../../reference'):
    """Make reference data into a given directory.

    :param dirname: Where to put the files.
    """
    try:
        os.mkdir(dirname)
    except OSError:
        pass

    print('Making Muller Data')
    make_muller_reference_data(dirname)

    print('Making Alanine Data')
    make_alanine_reference_data(dirname)


def parse():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate reference data for quant accel",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('dirname', help='Directory to write data',
                        default='./reference')

    args = parser.parse_args()
    make_reference_data(args.dirname)


if __name__ == "__main__":
    parse()