# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 14:19:53 2014

@author: harrigan
"""

import matplotlib
matplotlib.use('Agg')

import os
import numpy as np
from matplotlib import pyplot as pp
from matplotlib.colors import Normalize, LogNorm
import sys
import re
import scipy.io

from toy_accel import mullerforce as mf

from collections import defaultdict
import scipy.sparse.linalg
from scipy.sparse.linalg.eigen.arpack import ArpackNoConvergence
import logging as log

TEMP = 750
# Boltzmann constant in md units
KB = 0.0083145

class CentroidResult(object):
    def __init__(self):
        self.centroids_fn = None
        self.tmat_fn = None
        self.params = None
        self.round_i = None
        self.abspath = None

def _defaultdict_CentroidResult():
    return defaultdict(CentroidResult)

def walk(walkydir, centroid_regex, tmat_regex):
    # Two level defaltdict
    results = defaultdict(_defaultdict_CentroidResult)

    for dirpath, dirnames, filenames in os.walk(walkydir):
        log.info("Looking in directory %s", dirpath)
        for fn in filenames:
            # Regular expression to find files of interest
            centroid_match = re.match(centroid_regex, fn)
            tmat_match = re.match(tmat_regex, fn)
            complete_fn = os.path.join(dirpath, fn)

            # Whether there is any match
            match = None

            # Check if its a centroid file
            if centroid_match:
                match = centroid_match
                mtype = 'centroid'

            # Check if it's a tmat file
            if tmat_match:
                match = tmat_match
                mtype = 'tmat    '

            # If there is any match
            if match:
                # Get the round index
                round_i = int(match.group(1))

                # Hacky get params
                cwd = os.path.abspath(dirpath)
                param_str = cwd[cwd.rfind('/')+1:]
                fn_splits = param_str.split('-')
                try:
                    params = {fn_splits[0]: fn_splits[1]}
                except IndexError:
                    params = {param_str: param_str}

                # Debug statement
                log.debug("Found %s %d\tFilename: %s\tParams: %s",
                          mtype, round_i, complete_fn, params)

                # Get the result from the dictionary or make a new one
                result = results[param_str][round_i]
                if centroid_match:
                    # Load centroids with numpy
                    result.centroids_fn = complete_fn
                elif tmat_match:
                    # Load transition matrix with scipy
                    result.tmat_fn = complete_fn




                # In any event, save the params and round information
                result.params = params
                result.round_i = round_i
                result.abspath = cwd

    return results

def load(result):
    """Load the actual data given an object that has their filenames.
    """
    # Load centroids
    centroids = np.loadtxt(result.centroids_fn)

    # Load transition matrix
    tmat = scipy.io.mmread(result.tmat_fn)
    tmat = tmat.transpose()

    try:
        vals, vecs = scipy.sparse.linalg.eigs(tmat, k=1, which="LR",
                                              maxiter=100000, tol=1e-30)
        vecs = np.real_if_close(vecs)
        eq_distr = vecs[:, 0]
    except ArpackNoConvergence:
        log.warn("No eigenv convergence")
        eq_distr = np.ones(tmat.shape[0])

    # Compute a normalized equilibrium distribution
    eq_distr /= np.sum(eq_distr)

    return centroids, tmat, eq_distr


def project(frame, param_str):
    pass

def scatter(frame, param_str):
    """Scatter plot centroids where size and color are based on population.
    """
    # Make a scatter
    pp.subplot(121)
    pp.scatter(frame[0], frame[1], c=frame[2], s=2000*frame[2], norm=Normalize(vmin=0))
    pp.title(param_str)
    pp.colorbar()

    # Make theoretical distribution
    pp.subplot(122)
    calc_eq = mf.MullerForce.potential(frame[0], frame[1])
    calc_eq = np.exp(-calc_eq / (TEMP * KB))
    calc_eq /= np.sum(calc_eq)
    pp.scatter(frame[0], frame[1], c=calc_eq, s=2000*calc_eq, norm=Normalize(vmin=0))
    pp.title("Theoretical")
    pp.colorbar()


def make_movie(param_str, results, movie_dirname, movie='centroid'):
    """Make a movie for a given accelerator run.

        results: a dict of frames that has tmat_fn and centroids_fn for us
                 to load
        param_str: TODO
        movie_dirname: output folder name
        movie: ['centroid', 'projection']: Whether to make a movie from
               the centroids or project on to a grid
    """

    abs_movie_dirname = os.path.join(results.items()[0][1].abspath, movie_dirname % movie)
    log.info("Making %s movie for %s in directory %s", movie, param_str, abs_movie_dirname)

    # Make the directory
    try:
        os.mkdir(abs_movie_dirname)
    except OSError as e:
        log.warn(e)


    for round_i, frame_object in results.items():
        # Load up from file
        centroids, _, eq_distr = load(frame_object)
        if centroids is not None and eq_distr is not None:

            # Get relevant info into a list
            frame = [centroids[:, 0], centroids[:, 1], eq_distr]

            # Give a little debug info
            log.debug("Frame %d, centroid x's: %d, centroid y's: %d, eq_distr: %d",
                      round_i, len(frame[0]), len(frame[1]), len(frame[2]))

            # Plot
            if len(frame[0]) == len(frame[1]) and len(frame[0]) == len(frame[2]):
                if movie == 'centroid':
                    scatter(frame, param_str)
                elif movie == 'projection':
                    pass
                else:
                    log.error("Unknown movie type %s", movie)
            else:
                log.warn("Incompatible dimensions x: %d, y: %d, eq: %d",
                         len(frame[0]), len(frame[1]), len(frame[2]))
                pp.title("Incompatible.")
        else:
            # Something went wrong during loading
            if centroids is None:
                warning = "Centroids"
            elif eq_distr is None:
                warning = "Distribution"
            else:
                warning = "Something else"

            warning += " is none in frame %d for params %s"
            warning = warning % (round_i, param_str)
            pp.title(warning)
            log.warn(warning)

        # Save and clear
        fn_formatstr = os.path.join(abs_movie_dirname, '%s-%03d.png')

        pp.gcf().set_size_inches(14, 5)
        pp.savefig(fn_formatstr % ('frame', round_i))
        pp.clf()


def make_movies(all_results, movie_dirname):
    log.info("Making %d movies", len(all_results.keys()))
    log.debug("Different param configurations: %s", str(all_results.keys()))

    for param_str, subresults in all_results.items():
        make_movie(param_str, subresults, movie_dirname)


def main(argv):

    if len(argv) > 2 and argv[2] == 'percent':
        log.info("Using percent in regexp.")
        how = 'percent'
    else:
        how = 'round'

    results = walk(walkydir=argv[1],
                   centroid_regex='centroids-{how}-mkiii-([0-9]+).npy'.format(how=how),
                   tmat_regex='tmatfromclus-{how}-mkiii-([0-9]+).mtx'.format(how=how))

    make_movies(results, '%s-movie')

if __name__ == "__main__":
    log.basicConfig(level=log.INFO)
    main(sys.argv)
