# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:19:46 2014

@author: harrigan
"""

from quantaccel.tmat_simulation import RunResult
import pickle
import os
import re
import numpy as np
import scipy.optimize
from matplotlib import pyplot as pp
from collections import defaultdict

import logging as log


def load_runresults(rootdir):
    """From a directory, load all the results."""
    dirlist = os.listdir(rootdir)
    results = list()
    for fn in dirlist:
        if re.match("result-d-runcopy-[0-9]+-[0-9]+\\.pickl", fn):
            with open(os.path.join(rootdir, fn)) as f:
                r = pickle.load(f)
                results.append(fix_result(r))
    return results


def fix_result(r):
    """Process the results to make them fittable

    e.g. make real, remove outliers.
    """
    if not np.all(np.isreal(r.errors)):
        raise Exception("Unreal!")

    r.errors = np.delete(r.errors, np.where(r.errors[:, 1] > 1.0)[0], 0)
    return r


def ifunc(x, a):
    """Inverse."""
    return a / x


def ifunc_off(x, a, offset):
    """Inverse with a y-offset."""
    return a / x - offset


def simple_exp(x, tau):
    """Exponential with a y-offset."""
    return np.exp(-x / tau)

def bendy_exp(x, tau, n):
    """Exponential where the power of the exponent can vary."""
    return np.exp(- (x / tau) ** n)


def expifunc(x, a, n):
    """Inverse with tunable exponent."""
    return a / (x ** n)

def double_exp(x, c, d, tau1, tau2):
    return c * np.exp(-x/tau1) + d * np.exp(-x/tau2)


def add_result_to_dict(out_dict, speedup, quality, params):
    """After fitting, add it to some data structure.

    Note: no longer a dictionary.
    """
    to_append = (speedup, quality, params)
    out_dict.append(to_append)


def calc_speedup(fit_func, popt, err_bits=0.1, scale=1e6):
    """Calculate speedup for a given fit function and parameters

    err_bits is reference 'error'.
    """
    if fit_func == expifunc:
        a, n = popt
        return scale * (err_bits / a) ** (1.0 / n)
    elif fit_func == simple_exp:
        tau = popt
        return -scale * (1/(tau * np.log(err_bits)))
    elif fit_func == bendy_exp:
        tau, n = popt
        time = tau * (-np.log(err_bits))**(1.0/n)
        return scale / time
    else:
        raise Exception("No speedup calc for this fit func.")


def fit_results(results, fig_out_dir, fit_func=simple_exp, p0=None, show=False):
    """Fit all the results from a particular runcopy (directory of results)
    """
    n_fit = 0

    if p0 is None:
        if fit_func == expifunc:
            p0 = [1000., 1.]
        elif fit_func == simple_exp:
            p0 = [1000.]
        else:
            raise Exception("Please specify p0")

    # Sorry, it's no longer a dictionary
    out_dict = list()

    try:
        os.mkdir(fig_out_dir)
    except:
        pass

    for i, r in enumerate(results):
        pp.clf()
        pp.plot(r.errors[:, 0], r.errors[:, 1], 'o')
        try:
            if r.errors.shape[0] <= len(p0) * 2:
                raise Exception("Not enough data points")

            # Try to fit
            popt, pcov = scipy.optimize.curve_fit(
                fit_func, r.errors[:, 0], r.errors[:, 1], p0=p0)



            # Make fit curve
            if fit_func in [expifunc]:
                # For those that blow up at 0
                x = np.linspace(np.min(r.errors[:, 0]), np.max(r.errors[:, 0]))
            else:
                x = np.linspace(0, np.max(r.errors[:,0]))
            yfit = fit_func(x, *popt)
            pp.plot(x, yfit, '-')
            pp.title((str(r.params)))
            pp.xlabel(str(popt))

            # If the parameters do not change, it probably didn't fit
            if np.all(np.abs(popt - p0) < 1e-6):
                raise Exception("Same as p0")

            # We don't want any negative params either
            if np.any(popt<0):
                raise Exception("Negative params")

            speedup = calc_speedup(fit_func, popt)
            pp.suptitle("Speedup: %d" % int(speedup))
            # Add it to the dictionary
            add_result_to_dict(
                out_dict, speedup,
                1.0, r.params)
            n_fit += 1
        except Exception as e:
            pp.title(str(r.params))
            pp.suptitle("Didn't fit: %s" % e.args)

            add_result_to_dict(out_dict, 0, 0.0, r.params)

        pp.xlim(xmin=0)
        pp.ylim(ymin=0)
        pp.savefig(os.path.join(fig_out_dir, 'fit-%d.png' % i))
        if show: pp.show()
        log.debug("Fit result %d", i)

    log.info("Found fits for %d out of %i", n_fit, len(results))

    return out_dict


def average_many(superrootdir='.', regexp='d-run-[0-9]+'):
    """Fit over a set of runcopys (directories)."""

    # Get all our directories
    dirlist = os.listdir(superrootdir)
    dirlist = [fn for fn in dirlist if re.match(regexp, fn)]

    out_list = []

    for fn in dirlist:
        fullpath = os.path.join(superrootdir, fn)
        log.info("Loading results from %s", fullpath)
        results = load_runresults(rootdir=fullpath)
        out_list += fit_results(results, os.path.join(fullpath, 'figs/'))
    return out_list


def main():
    """call `average_many` and save the results."""
    out_dict = average_many()
    with open('d-run.list.pickl', 'w') as f:
        pickle.dump(out_dict, f)

if __name__ == "__main__":
    log.basicConfig(level=log.INFO)
    main()
