'''
Created on Mar 20, 2014

@author: harrigan
'''
from __future__ import division
import numpy as np
from quantaccel.tmat_simulation import RunResult
from quantaccel import fit_and_lump
import os
import pickle
import re
import collections
import scipy
import scipy.stats
import logging

log = logging.getLogger()
log.setLevel(logging.INFO)

NPOINTS = 50.0


def find_first_acceptable(r, bigger, cutoff):
    """Find the first point that is above or below a cutoff."""

    # Sort via time
    r.errors = r.errors[np.lexsort((r.errors[:, 1], r.errors[:, 0]))]

    # Find indices that match criterion
    if bigger:
        winds = np.where(r.errors[:, 1] > cutoff)[0]
    else:
        winds = np.where(r.errors[:, 1] < cutoff)[0]

    # pick the first
    if len(winds) > 0:
        ftime = r.errors[winds[0], 0]
        wind = winds[0]
    else:
        ftime = np.Inf
        wind, _ = r.errors.shape

    return ftime, wind


def avg_and_errorbar(values):
    """Take average in log space."""
    x_vals = sorted(set(values[:, 0]))
    avg_vals = np.zeros((len(x_vals), 4))

    for i, xval in enumerate(x_vals):
        selex = np.where(values[:, 0] == xval)[0]
        avg_vals[i, 0] = xval
        avg_vals[i, 1] = np.exp(np.mean(np.log(values[selex, 1])))

        log_std = np.std(np.log(values[selex, 1]))
        t_up = avg_vals[i, 1] * (np.exp(log_std) - 1)
        t_lo = avg_vals[i, 1] * -1 * (np.exp(-log_std) - 1)
        avg_vals[i, 3] = t_up
        avg_vals[i, 2] = t_lo
    return avg_vals


def _linear(x, m, b):
    """y = -mx + b."""
    return -m * x + b


def _exp(x, tau, c, off):
    """y = c * exp[-x/tau]"""
    return np.exp(-x / tau) * c + off


def _quad(x, a, b, c):
    """y = ax^2 + bx + c"""
    return a * x * x + b * x + c


def histogram_np(xval, y, hrange, n_bins, norm=False):
    """Use np.histogram to yield a normal histogram
    """
    hist, bin_edges = np.histogram(y, bins=n_bins, range=hrange, density=norm)
    bin_centers = np.diff(bin_edges) / 2.0 + bin_edges[:-1]
    return (hist, bin_centers, xval)


def histogram_ints(xval, y, hrange, n_bins):
    """Use np.bincount to yield counts for each integer.

    This works because number of rounds is an ostensibly an integer quantity
    """
    hist = np.bincount(y)
    bin_centers = np.arange(np.amax(y) + 1)
    return (hist, bin_centers, xval)


def histogram_kde(xval, y, hrange):
    """Use a kernel density estimator."""

    if hrange is not None:
        xmin, xmax = hrange
    else:
        xmin, xmax = (np.min(y), np.max(y))
    # Note! xvals is the output xaxis. These variable names are not good
    xvals = np.linspace(xmin, xmax, 100)

    if np.abs(np.max(y) - np.min(y)) < 1.0e-8:
        log.warn("%f is too skinny. Using counts", xval)
        return histogram_np(xval, y, hrange, n_bins=20, norm=True)

    kernel = scipy.stats.gaussian_kde(y)
    yvals = kernel.evaluate(xvals)
    return (yvals, xvals, xval)


def histogram(values, n_bins=15, whole_range=True, log=True, htype='np'):
    """Yield histograms for each 'x' """

    # Optionally take log of data
    if log:
        ally = np.log(values[:, 1])
    else:
        ally = np.array(values[:, 1], dtype=int)
    allx = values[:, 0]

    # Set range
    if whole_range:
        hrange = ((np.min(ally) - 1, np.max(ally) + 1))
    else:
        hrange = None

    # Pick function
    if htype == 'np':
        def hfunc(xval, y):
            return histogram_np(xval, y, hrange, n_bins)
    elif htype in ['int', 'ints']:
        def hfunc(xval, y):
            return histogram_ints(xval, y, hrange, n_bins)
    elif htype in ['kde']:
        def hfunc(xval, y):
            return histogram_kde(xval, y, hrange)
    else:
        raise ValueError('Specify a valid histogram type')

    # Do the main loop
    x_vals = sorted(set(allx))
    for i, xval in enumerate(x_vals):
        selex = np.where(allx == xval)[0]
        y = ally[selex]
        yield hfunc(xval, y)


class Results(object):
    """An object with all the results for a run."""

    def __init__(self):
        self.popresults = None
        self.itresults = None
        self.subdivide = None

        self.all_results = None
        self.na_results = None

        self.popresults_na = None

    def speed_pop(self, tvd_cutoff=0.6, name='pop'):
        """Find speed for population convergence, save results into the object.
        """
        for r in self._get_by_name(name):
            r.poptime, r.popind = find_first_acceptable(r, False, tvd_cutoff)

    def speed_subdivide_pop(self, tvd_cutoff=0.6, overwrite=False):
        """Take results that result from subdividing round 0 of those that
        converged before an adaptive round

        :param tvd_cutoff:
        :param overwrite: Whether to overwrite the index with the fractional one
        :return:
        """
        for r in self._get_by_name('sub'):
            r.poptime, r.popind = find_first_acceptable(r, False, tvd_cutoff)

        # Iterate through and match up in naive double-for loop
        for r in self._get_by_name('pop'):
            if r.popind == 0:
                for s in self.subdivide:
                    if r.params == s.params:
                        # Subtract 1 because later we add it on
                        if s.popind == 0:
                            log.warn("Subdivide converged in 0: %s",
                                     r.params)
                            s.popind = 0.1
                        r.spopind = (s.popind / s.errors.shape[0]) - 1
                        if overwrite:
                            r.popind = r.spopind
                        break
                else:
                    log.warn("Couldn't find subdivide: %s", r.params)

    def speed_it(self, it_cutoff=16605 / 1.5):
        """Find speed for it convergence, save results into the object."""
        for r in self.itresults:
            r.ittime, r.ittind = find_first_acceptable(r, True, it_cutoff)

    def get_runcopy(self, runcopy, name):
        """Select all the points for a specified runcopy."""
        sel_from = self._get_by_name(name)
        return [res for res in sel_from if res.params['runcopy'] == runcopy]

    def get_unique(self, name, *param_names):
        """Get unique param configurations."""

        # Apply compatibility
        param_names = [self.param_compat[pn] for pn in param_names]

        # Formatting string
        fstring = ["%s-%%s" % pn for pn in param_names if pn is not None]
        fstring = '_'.join(fstring)

        # Tuple generator
        def tuple_gen(r):
            return tuple([r.params[pn] for pn in param_names if pn is not None])

        # Set
        param_strs = set()
        for r in self._get_by_name(name):
            pvals = tuple_gen(r)
            pstring = fstring % pvals
            param_strs.add((pvals, pstring))

        param_strs = sorted(param_strs, key=lambda x: x[0])

        return fstring, tuple_gen, param_strs


    def plot_vs(self, name, yaxis, xaxis, label, where=None):
        """Get data in an (x,y) format for plotting.

        :name: type of convergence (pop, it)
        :yaxis: thing to show on yaxis (rounds, yval, [speedup])
        :xaxis: thing to show on xaxis (tpr, spt)
        :label: thing to use as labels

        ::returns:
            Dictionary of ndarrays shape (n, 2) of (xy) plot points. Key
            is the label
        """

        atlabel = PlotVS(name, yaxis, xaxis, label)

        if where is None:
            def where(r):
                return True

        for r in self._get_by_name(name):
            if yaxis == 'speedup':
                time = r.poptime
                lts = self._get_by_name(name, lts=True)
                yval = lts[r.params['tpr']] / time
            elif yaxis == 'yval':
                time = r.poptime
                yval = 1000.0 / time
            elif yaxis == 'rounds':
                yval = r.popind + 1
            elif yaxis == 'time':
                yval = r.poptime
            else:
                print "Not supported"

            if where(r):
                atlabel[r.params[label]] = np.append(
                    atlabel[r.params[label]], [[int(r.params[xaxis]), yval]],
                    axis=0)

        return atlabel


class MullerResults(Results):
    """A subclass for muller potential results."""

    def _get_by_name(self, name, lts=False):
        if name == 'pop' or name == 'projection-pop':
            if lts:
                return self.poplts
            else:
                return self.popresults
        elif name == 'it' or name == 'centroid-it':
            if lts:
                return self.itlts
            else:
                return self.itresults
        elif name == 'sub' or name == 'subdivide':
            return self.subdivide

    def load(self, file_list_fn, base_dir='.'):
        """Load results whose filename is listed in file_list_fn).

        :base_dir: filenames relative to this
        """
        results = list()
        with open(os.path.join(base_dir, file_list_fn)) as f:
            for line in f:
                line = line.strip()
                if len(line) == 0: continue
                with open(os.path.join(base_dir, line)) as pickl_f:
                    result = pickle.load(pickl_f)
                    match = re.search(r"runcopy-([0-9]+)", line)

                    # Try to figure out the runcopy number
                    if match is not None:
                        result.params['runcopy'] = int(match.group(1))

                    # Make everything into ints
                    for par in result.params:
                        try:
                            result.params[par] = int(result.params[par])
                        except ValueError:
                            pass
                    results.append(result)
        return results


    @property
    def param_compat(self):
        return {
            'lt': 'lt',
            'tpr': 'tpr',
            'spt': 'spt',
            'adaptive': None
        }


class TmatResults(Results):
    """A subclass for tmat results."""

    def _get_by_name(self, name, lts=False):

        # TODO: Take care of the fact that we don't use lts
        if lts:
            raise ValueError('No')

        if name == 'pop' or name == 'projection-pop':
            if self.popresults is None:
                self.popresults = list()
                for r in self.all_results:
                    rr = RunResult(r.params)
                    rr.errors = r.poperrors
                    self.popresults += [rr]
            return self.popresults
        if name == 'non' or name == 'non-adaptive':
            if self.popresults_na is None:
                self.popresults_na = list()
                for r in self.na_results:
                    rr = RunResult(r.params)
                    rr.errors = r.poperrors
                    self.popresults_na += [rr]
            return self.popresults_na


    def load(self, base_dir, letter):
        """Load results from a particular directory configuration.

        :letter: results have a file name with a specific letter.
        """

        run_dir_regex = r'%s-run-[0-9]+' % letter

        # list dirs e.g. k-run-21
        run_dirs = os.listdir(base_dir)
        run_dirs = [rd for rd in run_dirs if re.match(run_dir_regex, rd)]
        run_dirs = sorted(run_dirs,
                          key=lambda x: int(re.search('[0-9]+', x).group(0)))

        results = list()
        na_results = list()
        for rd in run_dirs:
            results += fit_and_lump.load_runresults(os.path.join(base_dir, rd),
                                                    adaptive=True)
            na_results += fit_and_lump.load_runresults(
                os.path.join(base_dir, rd), adaptive=False)

        self.all_results = results
        self.na_results = na_results
        log.info("Loaded %d points", len(results))

    @property
    def param_compat(self):
        return {
            'lt': None,
            'tpr': 'n_tpr',
            'spt': 'n_spt',
            'adaptive': 'adaptive'
        }


PRETTY_NAMES = {
    'spt': 'Trajectory length',
    'tpr': 'Parallelization',
    'lt': 'Lagtime',
    'adaptive': 'Adaptive',
    'rounds': '# Adaptive Rounds',
    'time': 'Time',
    None: 'Missing label'
}

PRETTY_NAMES['n_spt'] = PRETTY_NAMES['spt']
PRETTY_NAMES['n_tpr'] = PRETTY_NAMES['tpr']


class PlotVS(collections.defaultdict):
    """Subclass of dictionary to contain many plot lines where
    the label of the line is the key in the dictionary.

    Extended to support saving fits alongside, etc.
    """

    def __init__(self, name, yaxis, xaxis, label):
        super(PlotVS, self).__init__(lambda: np.zeros((0, 2)))
        self.fit_results = dict()
        self.xaxis = xaxis

        self.xlabel = PRETTY_NAMES[xaxis]
        self.ylabel = PRETTY_NAMES[yaxis]
        self.labellabel = PRETTY_NAMES[label]

    def items(self):
        """Return items sorted by key."""
        items = super(PlotVS, self).items()
        items = sorted(items, key=lambda x: int(x[0]))
        return items

    def enum_fit(self, norm_one=None, norm_all=None):
        """Yield a fit object for each possible fit.

        We defer performing the fit so plotting can happen."""

        if self.xaxis in ['spt', 'n_spt']:
            return self.enum_fit_vs_spt()
        elif self.xaxis in ['tpr', 'n_tpr']:
            return self.enum_fit_vs_tpr()

    def enum_fit_vs_spt(self):
        def dat_time_calc(x, y):
            return x * y

        for key, vals in self.items():
            fit = MullerFitResult(key, vals, dat_time_calc)
            yield fit

    def enum_fit_vs_tpr(self):
        for key, vals in self.items():
            fit = MullerFitResult(key, vals, lambda x, y: key * y)
            yield fit


class FitResult(object):
    """
    :param dat_time_calc: Take in x, y arrays and return array of time

    """

    def __init__(self, key, vals, dat_time_calc):
        """An object for a fit."""
        self.popt = None
        self.pcov = None
        self.x = vals[:, 0]
        self.y = vals[:, 1]
        self.fitx = None
        self.fity = None
        self.fit_time = None
        self.dat_time = dat_time_calc(self.x, self.y)
        self.key = key


class TmatFitResult(FitResult):
    pass


class MullerFitResult(FitResult):
    pass
