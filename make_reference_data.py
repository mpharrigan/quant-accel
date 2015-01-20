"""A script to generate reference data for testing.

This can be run from the command line.
"""

import os
from os.path import join as pjoin
import argparse
import pickle
import urllib.request
import tarfile
import shutil

import mdtraj.io
from mixtape.datasets.alanine_dipeptide import fetch_alanine_dipeptide
from maccelerator.configurations.muller import generate_muller_sysint
from maccelerator.configurations.alanine import generate_alanine_msm
from maccelerator.configurations.srckinase import generate_srckinase_msm
from maccelerator.simulate import OpenMMSimulator


SRC_URL = "https://stacks.stanford.edu/file/druid:cm993jk8755/"
SRC_FILE = "MSM_2000states_csrc.tar.gz"
SRC_DIR = "srckinase"


def get_src_kinase_data(dirname, cleanup=True):
    """Get the 2000 state msm from Stanford's SDR."""
    try:
        os.mkdir(pjoin(dirname, SRC_DIR))
    except OSError:
        pass

    # Fetch data
    with urllib.request.urlopen(SRC_URL + SRC_FILE) as tmat_tar_url:
        with open(pjoin(dirname, SRC_DIR, SRC_FILE), 'wb') as tmat_tar_f:
            tmat_tar_f.write(tmat_tar_url.read())
        with tarfile.open(pjoin(dirname, SRC_DIR, SRC_FILE)) as tmat_tar_f:
            tmat_tar_f.extractall(pjoin(dirname, SRC_DIR))

    # Load and convert
    msm, centers = generate_srckinase_msm(
        tmat_fn=pjoin(dirname, SRC_DIR, 'Data_l5', 'tProb.mtx'),
        populations_fn=pjoin(dirname, SRC_DIR, 'Data_l5', 'Populations.dat'),
        mapping_fn=pjoin(dirname, SRC_DIR, 'Data_l5', 'Mapping.dat'),
        gens_fn=pjoin(dirname, SRC_DIR, 'Gens.lh5')
    )

    # Save centers
    mdtraj.io.saveh(pjoin(dirname, 'src.centers.h5'),
                    cluster_centers=centers)

    # Save MSM Object
    with open(pjoin(dirname, 'src.msm.pickl'), 'wb') as f:
        pickle.dump(msm, f)

    # Optionally, delete all data
    if cleanup:
        shutil.rmtree(pjoin(dirname, SRC_DIR))


def make_alanine_reference_data(dirname):
    """Make a small transition matrix from Alanine trajectories

    We featurize using phi / psi angles

    Note: This function is not-deterministic, although it would be useful
    if it were, so testing could be conducted.
    """

    ala = fetch_alanine_dipeptide()
    msm, kmeans = generate_alanine_msm(ala)

    # Save cluster centers
    mdtraj.io.saveh(pjoin(dirname, 'ala.centers.h5'),
                    cluster_centers=kmeans.cluster_centers_)

    # Save MSM Object
    with open(pjoin(dirname, 'ala.msm.pickl'), 'wb') as f:
        pickle.dump(msm, f)


def make_muller_reference_data(dirname):
    """Make OpenMM system.xml and integrator.xml files.

    :param dirname: Where to put the files.
    """

    # TODO: Integrator has a random seed
    system, integrator = generate_muller_sysint()
    OpenMMSimulator.serialize(
        system, integrator,
        pjoin(dirname, 'muller_sys.xml'),
        pjoin(dirname, 'muller_int.xml')
    )


def make_reference_data(dirname, alanine=True, muller=True, srckinase=True):
    """Make reference data into a given directory.

    :param dirname: Where to put the files.
    """
    try:
        os.mkdir(dirname)
    except OSError:
        pass

    if alanine:
        print('Making Muller Data')
        make_muller_reference_data(dirname)

    if muller:
        print('Making Alanine Data')
        make_alanine_reference_data(dirname)

    if srckinase:
        print("Getting Src Kinase Data")
        get_src_kinase_data(dirname)


def parse():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate reference data for quant accel",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('dirname', help='Directory to write data',
                        default='./reference')

    # Flags
    parser.add_argument('--alanine', action='store_true', default=False)
    parser.add_argument('--muller', action='store_true', default=False)
    parser.add_argument('--srckinase', action='store_true', default=False)

    args = parser.parse_args()

    if (not args.alanine) and (not args.muller) and (not args.srckinase):
        make_reference_data(args.dirname)
    else:
        make_reference_data(args.dirname, alanine=args.alanine,
                            muller=args.muller, srckinase=args.srckinase)


if __name__ == "__main__":
    parse()
