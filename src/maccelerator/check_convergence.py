'''
Created on Mar 24, 2014

@author: harrigan

From an already built msm, check convergence and write to a file
'''

import argparse
from quantaccel import centroid_movie as cm
import os

TMAT_CENTROIDS = "../../spring_layout_centroids.npy"
PRECOMPUTED_EIGEN = "../../calc_eq.npy"

class ConvergenceChecker(object):
    pass

class PopulationProjectionTVD(ConvergenceChecker):
    pass


class PopulationCentroidTVD(ConvergenceChecker):
    pass



def main(args):
    """Take argparse args and call the function."""
    centroids_fn = 'centroids-{how}-mk{version}-{round_i}.npy'.format(
        **vars(args))
    tmat_fn = 'tmatfromclus-{how}-mk{version}-{round_i}.mtx'.format(
        **vars(args))
    mapping_fn = 'mapping-{how}-mk{version}-{round_i}.npy'.format(**vars(args))

    if args.system_type == 'muller':
        check_convergence_muller(args.round_i, centroids_fn, tmat_fn, 'msms')
    elif args.system_type == 'tmat':
        check_convergence_tmat(args.round_i, TMAT_CENTROIDS, tmat_fn,
                               mapping_fn, 'msms')


def write_convergence(round_i):
    """Write a file that indicates convergence was achieved."""
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


def check_convergence_tmat(round_i, centroids_fn, tmat_fn, mapping_fn,
                           dir_path):
    # Construct this object by hand
    result_obj = cm.CentroidResult()
    result_obj.centroids_fn = os.path.abspath(centroids_fn)
    result_obj.tmat_fn = os.path.abspath(os.path.join(dir_path, tmat_fn))
    result_obj.round_i = round_i
    result_obj.abspath = dir_path
    result_obj.mapping_fn = os.path.abspath(os.path.join(dir_path, mapping_fn))

    def helper(cen, tm, mapp):
        """Quick function to pass in the filename for the calculated
        equilibrium populations."""
        return cm.load_centroid_eqdistr_tmat(cen, tm, mapp,
                                             precomputed_eigen_fn=PRECOMPUTED_EIGEN)

    _, _, _, _, _, errorval = cm.load(result_obj,
                                      helper=helper)


    # TODO: Increase tolerance for ntl9
    # Use TVD of 0.4 for tmat
    if errorval < 0.2:
        write_convergence(round_i)


