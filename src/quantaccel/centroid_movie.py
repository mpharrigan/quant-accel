# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 14:19:53 2014

@author: harrigan
"""
from __future__ import division

import matplotlib
matplotlib.use('Agg')

import os
import numpy as np
from matplotlib import pyplot as pp
from matplotlib.colors import Normalize, LogNorm
import re
import scipy.io
import argparse
import pickle

from toy_accel import mullerforce as mf
from quantaccel.tmat_simulation import RunResult

from collections import defaultdict
import scipy.sparse.linalg
from scipy.sparse.linalg.eigen.arpack import ArpackNoConvergence
import scipy.interpolate
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
    with open(result.tmat_fn) as tmat_f:
        tmat = scipy.io.mmread(tmat_f)
        tmat = tmat.transpose()

        # Parse number of wallsteps
        tmat_f.seek(0)   # go to comment line
        tmat_f.readline()  # read it
        try:
            # Parse the comment line
            wallsteps = int(tmat_f.readline().strip().split()[-1])
        except ValueError:
            wallsteps = -1
    log.debug("Number of wallsteps: %d", wallsteps)

    # Get the eigenvectors
    try:
        vals, vecs = scipy.sparse.linalg.eigs(tmat, k=1, which="LR",
                                              maxiter=100000, tol=1e-30)
        vecs = np.real_if_close(vecs)

        if np.abs(np.real(vals)-1) > 1e-8:
            raise ValueError('Eigenvalue is not 1')

        eq_distr = vecs[:, 0]
    except (ArpackNoConvergence, ValueError) as e:
        log.warn("No eigenv convergence %s", str(e))
        eq_distr = np.ones(tmat.shape[0])

    # Compute a normalized equilibrium distribution
    eq_distr /= np.sum(eq_distr)

    return centroids, tmat, eq_distr, wallsteps


def get_grid(volume, resolution):
    """Grid up a space

    resolution - how fine the grid should be
    """
    (xmin, xmax, ymin, ymax) = volume.bounds
    log.debug('Making grid with bounds %.2f %.2f %.2f %.2f', *volume.bounds)
    grid_width = max(xmax - xmin, ymax - ymin) / resolution
    log.debug("Gridwithd %f", grid_width)
    grid = np.mgrid[xmin : xmax : grid_width, ymin : ymax : grid_width]
    return grid


def distribution_norm_tvd(p, q):
    """Total variation distance.
    """

    q /= np.sum(q)
    p /= np.sum(p)

    res = 0.5 * np.sum(np.abs(p - q))

    return res

def distribution_norm(p, q):
    """Take the norm of a distribution.

        p: actual
        q: estimated
    """
    return distribution_norm_tvd(p, q)

def project(frame, param_str, grid):
    xx, yy = grid
    bounds = (xx.min(), xx.max(), yy.min(), yy.max())

    # Make the projection
    pp.subplot(121)
    pp.title(param_str)
    known_points = np.vstack((frame[0], frame[1])).T
    project_on = np.vstack([xx.ravel(), yy.ravel()]).T # todo: remove
    est = scipy.interpolate.griddata(known_points, frame[2],
                                     (xx, yy), method='cubic',
                                     fill_value = 0.0)
    est = est.clip(min=0.0)
    est /= np.sum(est)

    #est_show = (-np.log(est)).clip(max=70)
    est_show = est
    pp.imshow(est_show.T, interpolation='nearest',
              extent=bounds,
              aspect='auto',
              origin='lower')
    pp.colorbar()

    pp.subplot(122)
    pp.title("Theoretical")
    calc_eq = mf.MullerForce.potential(xx, yy)
    calc_eq = np.exp(-calc_eq / (TEMP * KB))
    calc_eq /= np.sum(calc_eq)
    #calc_eq_show = (-np.log(calc_eq)).clip(max=70)
    calc_eq_show = calc_eq
    pp.imshow(calc_eq_show.T, interpolation='nearest',
              extent=bounds,
              aspect='auto',
              origin='lower')
    pp.colorbar()

    return distribution_norm(calc_eq, est)

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

class Volume(object):
    def __init__(self, centroidx=None, centroidy=None):
        if centroidx is not None and centroidy is not None:
            self.xmin = np.min(centroidx)
            self.xmax = np.max(centroidx)
            self.ymin = np.min(centroidy)
            self.ymax = np.max(centroidy)
        else:
            self.xmin = 0.0
            self.xmax = 0.0
            self.ymin = 0.3
            self.ymax = 0.3

    @property
    def volume(self):
        return (self.xmax - self.xmin) * (self.ymax - self.ymin)

    @property
    def bounds(self):
        return (self.xmin, self.xmax, self.ymin, self.ymax)

    def union(self, other_v):
        """Compare to other volume and turn me into something that contains
        both.
        """
        self.xmin = min(self.xmin, other_v.xmin)
        self.ymin = min(self.ymin, other_v.ymin)
        self.xmax = max(self.xmax, other_v.xmax)
        self.ymax = max(self.ymax, other_v.ymax)

VOLUMES = dict()
BIGV = Volume()
SETVOL = Volume([-3.0, 1.0], [-1.0, 3.0])

def make_movie(param_str, results, movie_dirname, movie):
    """Make a movie for a given accelerator run.

        results: a dict of frames that has tmat_fn and centroids_fn for us
                 to load
        param_str: This goes in the title
        movie_dirname: output folder name
        movie: ['centroid', 'projection']: Whether to make a movie from
               the centroids or project on to a grid
    """
    assert movie in ['centroid', 'projection'], 'Invalid movie type'

    abs_movie_dirname = os.path.join(results.items()[0][1].abspath, movie_dirname % movie)
    log.info("Making %s movie for %s in directory %s", movie, param_str, abs_movie_dirname)

    # Make the directory
    try:
        os.mkdir(abs_movie_dirname)
    except OSError as e:
        log.warn(e)

    biggest_v = Volume()

    errors = -1 * np.ones((len(results.keys()), 2))

    if movie == 'projection':
        grid = get_grid(SETVOL, resolution=200)

    for round_i, frame_object in results.items():
        # Load up from file
        centroids, _, eq_distr, walltime = load(frame_object)
        if centroids is not None and eq_distr is not None:

            # Get relevant info into a list
            frame = [centroids[:, 0], centroids[:, 1], eq_distr]
            biggest_v.union(Volume(frame[0], frame[1]))

            # Give a little debug info
            log.debug("Frame %d, centroid x's: %d, centroid y's: %d, eq_distr: %d",
                      round_i, len(frame[0]), len(frame[1]), len(frame[2]))

            # Plot
            if len(frame[0]) == len(frame[1]) and len(frame[0]) == len(frame[2]):
                if movie == 'centroid':
                    scatter(frame, param_str)
                elif movie == 'projection':
                    dist = project(frame, param_str, grid)
                    errors[round_i - 1] = [walltime, dist]
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

    VOLUMES[param_str] = biggest_v
    BIGV.union(biggest_v)

    # Delete anything with -1 (no round)
    errors = np.delete(errors, np.where(errors[:,0] < 0)[0], axis=0)

    with open('quant-%s-%s.pickl' % (movie, param_str), 'w') as f:
        (_, someresult) = results.items()[0]
        pickle.dump(RunResult(someresult.params, errors), f, protocol=2)

    log.info("The largest volume for that movie is ((%.2f, %.2f), (%.2f, %.2f)) -> %.3f",
             *(biggest_v.bounds + (biggest_v.volume,)))


def make_movies(all_results, movie_dirname, movie_type):
    log.info("Making %d movies", len(all_results.keys()))
    log.debug("Different param configurations: %s", str(all_results.keys()))

    for param_str, subresults in all_results.items():
        make_movie(param_str, subresults, movie_dirname, movie_type)

    log.info("The largest volume overall is ((%.2f, %.2f), (%.2f, %.2f)) -> %.3f",
             *(BIGV.bounds + (BIGV.volume,)))


def parse():
    parser = argparse.ArgumentParser(description="Make movies etc")
    parser.add_argument('walkydir',
                        help='''initial dir to start walk from''')
    parser.add_argument('version',
                        help='''version. Will look for files with mk{version}''',
                        type=int)
    parser.add_argument('-how', dest='how',
                        help='''either percent or round''',
                        default='round')
    parser.add_argument('-mt', dest='movietype',
                        help='''Type of movie: [centroid, projections]''',
                        default='centroid')
    parser.add_argument('--debug', dest='debug',
                        help='''Print extremely verbose output''',
                        action='store_true')
    parser.set_defaults(debug=False)

    args = parser.parse_args()
    if args.debug: log.basicConfig(level=log.DEBUG)
    else: log.basicConfig(level=log.INFO)
    main(args.walkydir, args.how, args.version, args.movietype)

def main(walkydir, how, version, movietype):

    fmt = {'how': how, 'version': version}
    centroid_regex='centroids-{how}-mk{version}-([0-9]+).npy'.format(**fmt)
    tmat_regex='tmatfromclus-{how}-mk{version}-([0-9]+).mtx'.format(**fmt)

    results = walk(walkydir=walkydir,
                   centroid_regex=centroid_regex,
                   tmat_regex=tmat_regex)

    make_movies(results, '%s-movie-mk{version}'.format(**fmt), movietype)

if __name__ == "__main__":
    parse()
