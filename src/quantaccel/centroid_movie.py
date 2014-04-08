# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 14:19:53 2014

@author: harrigan
"""


from __future__ import division

import logging as log
log.basicConfig(level=log.INFO)

import matplotlib
matplotlib.use('Agg')

from collections import defaultdict
from matplotlib import pyplot as pp
from matplotlib.colors import Normalize, LogNorm
from msmbuilder import msm_analysis as msma
from quantaccel.tmat_simulation import RunResult
from scipy.sparse.linalg.eigen.arpack import ArpackNoConvergence
from toy_accel import mullerforce as mf
import argparse
import numpy as np
import os
import pickle
import re
import scipy.interpolate
import scipy.io
import scipy.sparse.linalg

TEMP = 750
# Boltzmann constant in md units
KB = 0.0083145
LAGTIME = 20


class CentroidResult(object):

    """Container for filenames and params."""

    def __init__(self):
        self.centroids_fn = None
        self.tmat_fn = None
        self.params = None
        self.round_i = None
        self.abspath = None


def _defaultdict_CentroidResult():
    """For making a defaultdict of defaultdicts."""
    return defaultdict(CentroidResult)


def walk(walkydir, centroid_regex, tmat_regex):
    """Walk over the directory structure looking for centroid and tmat files.
    """
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
                dnames = cwd.split('/')
                dname = dnames[-2]
                param_strs = dname.split('_')
                params = dict()
                for param_str in param_strs:
                    splits = param_str.split('-')
                    try:
                        params[splits[0]] = splits[1]
                    except IndexError:
                        params[param_str] = param_str
                

                # Debug statement
                log.debug("Found %s %d\tFilename: %s\tParams: %s",
                          mtype, round_i, complete_fn, params)

                # Get the result from the dictionary or make a new one
                result = results[dname][round_i]
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


def load(result, helper):
    """Load the actual data given an object that has their filenames.
    """
    # Load centroids
    centroids = np.loadtxt(result.centroids_fn)
    if len(centroids.shape) == 1:
        centroids = centroids[np.newaxis, :]

    # Load transition matrix
    with open(result.tmat_fn) as tmat_f:
        tmat = scipy.io.mmread(tmat_f)

        # Parse number of wallsteps
        tmat_f.seek(0)  # go to comment line
        tmat_f.readline()  # read it
        try:
            # Parse the comment line
            wallsteps = int(tmat_f.readline().strip().split()[-1])
        except ValueError:
            wallsteps = -1
    log.debug("Number of wallsteps: %d", wallsteps)
    est_distr, theory_distr, errorval = helper(centroids, tmat)
    return centroids, tmat, wallsteps, est_distr, theory_distr, errorval


def load_numstates(centroids, tmat):
    """Just put the number of states in the errorval."""
    return None, None, len(centroids)


def load_centroid_eigen(centroids, tmat):
    """Plot first eigenvalue on centroids and return errorval = 1st IT."""
    # Get the eigenvectors
    try:
        n = tmat.shape[0]
        if n < 2:
            raise ValueError('Transition matrix too small')
        elif n <= 3:
            tmat = tmat.todense()
        vals, vecs = msma.get_eigenvectors(tmat, n_eigs=2)
        if np.abs(vals[0] - 1) > 1e-6:
            log.warn("First eigenvalue is not one!")
        errorval = -LAGTIME / np.log(vals[1])
    except (ArpackNoConvergence, ValueError) as e:
        log.warn("No eigenv convergence %s", str(e))
        vecs = np.ones((tmat.shape[0], 2))
        errorval = 0

    est_plot = [centroids[:, 0], centroids[:, 1], vecs[:, 1]]
    # TODO: Get actual theoretical eigenvec
    theory_plot = [centroids[:, 0], centroids[:, 1], np.ones(tmat.shape[0])]

    return est_plot, theory_plot, errorval


def load_project_eqdistr(centroids, tmat):
    """Compute equilibrium density."""
    xx, yy = GRID
    n_states = tmat.shape[0]

    # Get the eigenvectors
    tmat = tmat.transpose()
    try:
        vals, vecs = scipy.sparse.linalg.eigs(tmat, k=1, which="LR",
                                              maxiter=100000, tol=1e-30)
    except (ArpackNoConvergence, ValueError) as e:
        log.warn('Sparse solver threw error: %s', str(e))
        vals, vecs = scipy.linalg.eig(tmat.toarray())
        
    order = np.argsort(-np.real(vals))
    vals = np.real_if_close(vals[order])
    vecs = np.real_if_close(vecs[:, order])

    if np.abs(np.real(vals[0]) - 1) > 1e-8:
        log.warn('Eigenvalue is not 1: %f', vals[0])
        
    
    eq_distr = vecs[:, 0]
    # Compute a normalized equilibrium distribution
    eq_distr /= np.sum(eq_distr)
    
    if len(centroids.shape) == 1:
        centroids = centroids[np.newaxis, :]

    known_points = np.vstack((centroids[:, 0], centroids[:, 1])).T
    if n_states >= 4: method = 'cubic'
    else: method = 'nearest'
    est = scipy.interpolate.griddata(known_points, eq_distr,
                                     (xx, yy), method=method,
                                     fill_value=0.0)
    est = est.clip(min=0.0)
    est /= np.sum(est)

    calc_eq = mf.MullerForce.potential(xx, yy)
    calc_eq = np.exp(-calc_eq / (TEMP * KB))
    calc_eq /= np.sum(calc_eq)

    errorval = distribution_norm(calc_eq, est)

    return est, calc_eq, errorval


def load_centroid_eqdistr(centroids, tmat):
    """Compute equilibrium density for centroids."""

    # Get the eigenvectors
    tmat = tmat.transpose()
    try:
        vals, vecs = scipy.sparse.linalg.eigs(tmat, k=1, which="LR",
                                              maxiter=100000, tol=1e-30)
        vecs = np.real_if_close(vecs)

        if np.abs(np.real(vals) - 1) > 1e-8:
            raise ValueError('Eigenvalue is not 1')

        eq_distr = vecs[:, 0]
    except (ArpackNoConvergence, ValueError) as e:
        log.warn("No eigenv convergence %s", str(e))
        eq_distr = np.ones(tmat.shape[0])

    # Compute a normalized equilibrium distribution
    # eq_distr = eq_distr.clip(min=0.0)
    eq_distr /= np.sum(eq_distr)

    calc_eq = mf.MullerForce.potential(centroids[:, 0], centroids[:, 1])
    calc_eq = np.exp(-calc_eq / (TEMP * KB))
    calc_eq /= np.sum(calc_eq)

    errorval = distribution_norm(calc_eq, eq_distr)

    est_plot = [centroids[:, 0], centroids[:, 1], eq_distr]
    theory_plot = [centroids[:, 0], centroids[:, 1], calc_eq]

    return est_plot, theory_plot, errorval


def load_it_distr(centroids, tmat):
    """Compute IT error and plot first eigen."""
    # TODO: Implement
    pass


def get_grid(volume, resolution):
    """Grid up a space

    :resolution: how fine the grid should be
    """
    (xmin, xmax, ymin, ymax) = volume.bounds
    log.debug('Making grid with bounds %.2f %.2f %.2f %.2f', *volume.bounds)
    grid_width = max(xmax - xmin, ymax - ymin) / resolution
    log.debug("Gridwithd %f", grid_width)
    grid = np.mgrid[xmin: xmax: grid_width, ymin: ymax: grid_width]
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

        :p: actual
        :q: estimated
    """
    return distribution_norm_tvd(p, q)


def project(est_plot, theory_plot, param_str):
    """Plot a projection onto a grid."""
    grid = get_grid(SETVOL, resolution=200)
    xx, yy = grid
    bounds = (xx.min(), xx.max(), yy.min(), yy.max())

    # Make the projection
    pp.subplot(121)
    pp.title(param_str)
    pp.imshow(est_plot.T, interpolation='nearest',
              extent=bounds,
              aspect='auto',
              origin='lower')
    pp.colorbar()

    pp.subplot(122)
    pp.title("Theoretical")
    pp.imshow(theory_plot.T, interpolation='nearest',
              extent=bounds,
              aspect='auto',
              origin='lower')
    pp.colorbar()


def scatter(est_plot, theory_plot, param_str, do_size):
    """Scatter plot centroids where size and color are based on population.
    """
    if do_size:
        est_s = 2000 * est_plot[2]
        theory_s = 2000 * theory_plot[2]
        norm = Normalize(vmin=0)
    else:
        est_s = 100
        theory_s = 100
        extent = max(np.max(est_plot[2]), -np.min(est_plot[2]))
        norm = Normalize(vmin=-extent, vmax=extent)

    # import pdb; pdb.set_trace()

    # Make a scatter
    pp.subplot(121)
    pp.scatter(est_plot[0], est_plot[1],
               c=est_plot[2], s=est_s, norm=norm)
    pp.title(param_str)
    pp.colorbar()

    # Make theoretical distribution
    pp.subplot(122)
    pp.scatter(theory_plot[0], theory_plot[1],
               c=theory_plot[2], s=theory_s, norm=norm)
    pp.title("Theoretical")
    pp.colorbar()


class Volume(object):

    """An object representing a rectangular area

    Use this to keep track of how big of a grid we need to make to
    project onto without discarding any info.
    """

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
        """The volume."""
        return (self.xmax - self.xmin) * (self.ymax - self.ymin)

    @property
    def bounds(self):
        """Bounds of this volume."""
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
#SETVOL = Volume([-3.0, 1.0], [-1.0, 3.0])
#minx=-1.5, maxx=1.2, miny=-0.2, maxy=2
SETVOL = Volume([-1.5, 1.2], [-0.2, 3.0])
GRID = get_grid(SETVOL, resolution=200)


def make_movie(param_str, results, movie_dirname, movie):
    """Make a movie for a given accelerator run.

        :results: a dict of frames that has tmat_fn and centroids_fn for us
                 to load
        :param_str: This goes in the title
        :movie_dirname: output folder name
        :movie: ['centroid-pop', 'projection-pop',
            'centroid-it', 'projection-it']: Whether to make a movie from
            the centroids or project on to a grid
    """

    if movie == 'centroid-pop':
        load_helper = load_centroid_eqdistr
        scatter_helper = lambda e, t, p: scatter(e, t, p, do_size=True)
    elif movie == 'projection-pop':
        load_helper = load_project_eqdistr
    elif movie == 'centroid-it':
        load_helper = load_centroid_eigen
        scatter_helper = lambda e, t, p: scatter(e, t, p, do_size=False)
    elif movie == 'projection-it':
        load_helper = None  # TODO
    elif movie == 'num-states':
        load_helper = load_numstates
    else:
        log.error("Invalid movie type %s", movie)

    abs_movie_dirname = os.path.join(
        results.items()[0][1].abspath,
        movie_dirname %
        movie)
    log.info(
        "Making %s movie for %s in directory %s",
        movie,
        param_str,
        abs_movie_dirname)

    # Make the directory
    try:
        os.mkdir(abs_movie_dirname)
    except OSError as e:
        log.warn(e)

    biggest_v = Volume()
    errors = -1 * np.ones((len(results.keys()), 2))

    for round_i, frame_object in results.items():
        # Load up from file
        # (centroids, tmat, walltime, est_distr, theory_distr, errorval)
        centroids, _, walltime, est_distr, theory_distr, errorval = \
            load(frame_object, helper=load_helper)
        errors[round_i - 1] = [walltime, errorval]

        biggest_v.union(Volume(centroids[:, 0], centroids[:, 1]))

        save_fig = True
        if 'centroid' in movie:
            scatter_helper(est_distr, theory_distr, param_str)
        elif 'projection' in movie:
            project(est_distr, theory_distr, param_str)
        else:
            log.warn("Not making a movie for %s", movie)
            save_fig = False

        # Save and clear
        fn_formatstr = os.path.join(abs_movie_dirname, '%s-%03d.png')

        if save_fig:
            pp.gcf().set_size_inches(14, 5)
            pp.savefig(fn_formatstr % ('frame', round_i))
            pp.clf()

    VOLUMES[param_str] = biggest_v
    BIGV.union(biggest_v)

    # Delete anything with -1 (no round)
    errors = np.delete(errors, np.where(errors[:, 0] < 0)[0], axis=0)

    with open('quant-%s-%s.pickl' % (movie, param_str), 'w') as f:
        (_, someresult) = results.items()[0]
        rr = RunResult(someresult.params)
        rr.errors = errors
        pickle.dump(rr, f, protocol=2)

    log.info("The largest volume for that movie is ((%.2f, %.2f), (%.2f, %.2f)) -> %.3f",
             *(biggest_v.bounds + (biggest_v.volume,)))


def make_movies(all_results, movie_dirname, movie_type):
    """Iterate through all folders and results in folders,
    make movie in :movie_dirname: of type :movie_type:
    """
    log.info("Making %d movies", len(all_results.keys()))
    log.debug("Different param configurations: %s", str(all_results.keys()))

    for param_str, subresults in all_results.items():
        make_movie(param_str, subresults, movie_dirname, movie_type)

    log.info("The largest volume overall is ((%.2f, %.2f), (%.2f, %.2f)) -> %.3f",
             *(BIGV.bounds + (BIGV.volume,)))


def parse():
    """Parse command line options."""
    parser = argparse.ArgumentParser(description="Make movies etc",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('walkydir',
                        help='''initial dir to start walk from''')
    parser.add_argument('version',
                        help='''version. Will look for files with mk{version}''',
                        type=int)
    parser.add_argument('-how', dest='how',
                        help='''either percent or round''',
                        default='round')
    parser.add_argument('--debug', dest='debug',
                        help='''Print extremely verbose output''',
                        action='store_true')
    parser.set_defaults(debug=False)
    parser.add_argument('--movietype', '-mt', choices=
                        ['centroid-pop', 'projection-pop', 'centroid-it',
                         'projection-it', 'num-states'],
                        help='''Type of movie''',
                        default='projection-pop')

    args = parser.parse_args()
#     if args.debug:
#         log.basicConfig(level=log.DEBUG)
#     else:
#         log.basicConfig(level=log.INFO)
    main(args.walkydir, args.how, args.version, args.movietype)


def main(walkydir, how, version, movietype):
    """Entry point"""

    fmt = {'how': how, 'version': version}
    centroid_regex = 'centroids-{how}-mk{version}-([0-9]+).npy'.format(**fmt)
    tmat_regex = 'tmatfromclus-{how}-mk{version}-([0-9]+).mtx'.format(**fmt)

    results = walk(walkydir=walkydir,
                   centroid_regex=centroid_regex,
                   tmat_regex=tmat_regex)

    make_movies(results, '%s-mk{version}'.format(**fmt), movietype)

if __name__ == "__main__":
    log.info("Logging is working")
    parse()
