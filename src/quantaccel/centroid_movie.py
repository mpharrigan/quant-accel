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
import sys
import re
import scipy.io
import shutil

from collections import defaultdict
from operator import attrgetter

import scipy.sparse.linalg
from scipy.sparse.linalg.eigen.arpack import ArpackNoConvergence

import logging as log

class CentroidResult(object):
    def __init__(self):
        self.eq_distr = None
        self.centroids = None
        self.params = None
        self.round_i = None
        self.abspath = None

def _defaultdict_CentroidResult():
    return defaultdict(CentroidResult)

def load(walkydir, centroid_regex, tmat_regex):
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
                params = {fn_splits[0]: fn_splits[1]}

                # Debug statement
                log.debug("Found %s %d\tFilename: %s\tParams: %s",
                          mtype, round_i, complete_fn, params)

                # Get the result from the dictionary or make a new one
                result = results[param_str][round_i]
                if centroid_match:
                    # Load centroids with numpy
                    centroids = np.loadtxt(complete_fn)
                    result.centroids = centroids
                elif tmat_match:
                    # Load transition matrix with scipy
                    tmat = scipy.io.mmread(complete_fn)

                    if result.centroids is not None:
                        centroidlen = centroids.shape[0]
                        tmatlen = tmat.shape[0]
                        if centroidlen != tmatlen:
                            log.warn("Problem with %s matching lenghts %d and %d",
                                     complete_fn, tmatlen, centroidlen)

                    try:
                        vals, vecs = scipy.sparse.linalg.eigs(tmat, k=1, which="LR", maxiter=100000, tol=1e-30)
                        vecs = np.real_if_close(vecs)
                        eq_distr = vecs[:,0]
                    except ArpackNoConvergence:
                        log.warn("No eigenv convergence")
                        eq_distr = np.ones(tmat.shape[0])

                    # eq_distr /= np.sum(eq_distr)
                    result.eq_distr = eq_distr


                # In any event, save the params and round information
                result.params = params
                result.round_i = round_i
                result.abspath = cwd

    return results

def make_movie(param_str, results, movie_dirname):
    abs_movie_dirname = os.path.join(results.items()[0][1].abspath, movie_dirname)
    log.info("Making movie for %s in directory %s", param_str, abs_movie_dirname)

    try:
        os.mkdir(abs_movie_dirname)
    except Exception as e:
        log.warn(e)

    plottable = list()
    for round_i, frame_object in results.items():
        # For each frame make a plot
        if frame_object.centroids is not None and frame_object.eq_distr is not None:

            # Get relevant info
            frame = [frame_object.centroids[:,0],
                     frame_object.centroids[:,1],
                     frame_object.eq_distr]

            # Give a little debug info
            log.debug("Frame %d, centroid x's: %d, centroid y's: %d, eq_distr: %d",
                      round_i, len(frame[0]), len(frame[1]), len(frame[2]))

            # Plot
            if len(frame[0]) == len(frame[1]) and len(frame[0]) == len(frame[2]):
                pp.scatter(frame[0], frame[1], s=500*frame[2], c=frame[2])
                pp.title(param_str)
                pp.colorbar()
            else:
                log.warn("Incompatible dimensions x: %d, y: %d, eq: %d",
                         len(frame[0]), len(frame[1]), len(frame[2]))
                pp.title("Incompatible.")
        else:
            # Something went wrong during loading
            if frame_object.centroids is None:
                warning = "Centroids"
            elif frame_object.eq_distr is None:
                warning = "Distribution"
            else:
                warning = "Something else"

            warning += " is none in frame %d for params %s"
            warning = warning % (round_i, param_str)
            pp.title(warning)
            log.warn(warning)

        # Save and clear
        pp.savefig(os.path.join(abs_movie_dirname, 'frame-%03d.png' % round_i))
        pp.clf()





def make_movies(all_results, movie_dirname):
    log.info("Making %d movies", len(all_results.keys()))
    log.debug("Different param configurations: %s", str(all_results.keys()))

    for param_str, subresults in all_results.items():
        make_movie(param_str, subresults, movie_dirname)


def main(argv):
    log.basicConfig(level=log.INFO)

    results = load(walkydir=argv[1],
                   centroid_regex='centroids-round-mkii-([0-9]+).npy',
                   tmat_regex='tmatfromclus-round-mkii-([0-9]+).mtx')

    make_movies(results, movie_dirname='centroid-movie-mki/')

if __name__ == "__main__":
    main(sys.argv)
