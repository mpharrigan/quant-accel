'''
Created on Mar 24, 2014

@author: harrigan

From an already built msm, check convergence and write to a file
'''

import argparse
from quantaccel import centroid_movie as cm
import os


def main(args):
    """Take argparse args and call the function."""
    centroids_fn = 'centroids-{how}-mk{version}-{round_i}.npy'.format(
        **vars(args))
    tmat_fn = 'tmatfromclus-{how}-mk{version}-{round_i}.mtx'.format(
        **vars(args))

    if args.system_type == 'muller':
        check_convergence_muller(args.round_i, centroids_fn, tmat_fn, 'msms')
    elif args.system_type =='tmat':
        check_convergence_tmat(args.round_i, tmat_fn, 'msms')

def write_convergence(round_i):
    with open('converged', 'w') as f:
        f.write("Converged at round %d\n" % round_i)

def check_convergence_muller(round_i, centroids_fn, tmat_fn, dir_path):
    """For a given tmat and centroids, check convergence.

    This uses equilibrium distribution. Maybe later we can add
    an option to use a different criterion.
    """

    # Construct this object by hand so we can use the code in centroid_movie
    result_obj = cm.CentroidResult()
    result_obj.centroids_fn = os.path.abspath(
        os.path.join(dir_path, centroids_fn))
    result_obj.tmat_fn = os.path.abspath(os.path.join(dir_path, tmat_fn))
    result_obj.round_i = round_i
    result_obj.abspath = dir_path
    _, _, _, _, _, errorval = cm.load(result_obj,
                                      helper=cm.load_project_eqdistr)

    # Use TVD of 0.6 for muller
    if errorval < 0.6:
        write_convergence(round_i)


def check_convergence_tmat(round_i, tmat_fn, dir_path):
    #TODO: Make movie code for tmat
    pass


def parse():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Check the convergence of an msm.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('system_type', choices=['muller', 'tmat'])

    parser.add_argument('round_i', type=int,
                        help="""The round to check convergence.""")
    parser.add_argument('how',
                        choices=['rnew'])
    parser.add_argument('version',
                        help="""Look for mk{version}""", type=int)

    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    parse()
