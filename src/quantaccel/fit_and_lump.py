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
        if re.match("result-d-runcopy-[0-9]+-[0-9]+\.pickl", fn):
            with open(os.path.join(rootdir, fn)) as f:
                r = pickle.load(f)
                results.append(fix_result(r))
    return results

def fix_result(r):
    """Process the results to make them fittable

    e.g. make real, remove outliers.
    """

    r.errors = np.delete(r.errors, np.where(r.errors[:,1] > 1)[0], 0)
    return r

def ifunc(x,a):
    """Inverse."""
    return a/x
def ifunc_off(x, a, offset):
    """Inverse with a y-offset."""
    return a/x - offset
def exp(x, tau, c, offset):
    """Exponential with a y-offset."""
    return c*np.exp(-x/tau) + offset
def expifunc(x, a, n):
    """Inverse with tunable exponent."""
    return a/(x**n)

def add_result_to_dict(out_dict, speedup, quality, params):
    to_append = (speedup, quality, params)
    out_dict.append(to_append)

def speedup(fit_func, popt, err_bits=0.1):
    """Calculate speedup for a given fit function and parameters

    err_bits is reference 'error'.
    """
    if fit_func == expifunc:
        a, n = popt
        return 1e6 * (err_bits / a)**(1.0/n)

def fit_results(results, fig_out_dir, fit_func=expifunc, p0=[1000., 1.]):
    n_fit = 0

    # Sorry, it's no longer a dictionary
    out_dict = list()

    try:
        os.mkdir(fig_out_dir)
    except:
        pass

    for i, r in enumerate(results):
        pp.clf()
        pp.plot(r.errors[:,0], r.errors[:,1], 'o')
        try:
            # Try to fit
            popt, pcov = scipy.optimize.curve_fit(fit_func, r.errors[:,0], r.errors[:,1], p0=p0)

            # If the parameters do not change, it probably didn't fit
            if np.all(np.abs(popt-p0)<1e-6): raise Exception

            # Make fit curve
            x = np.linspace(np.min(r.errors[:,0]), np.max(r.errors[:,0]))
            yfit = fit_func(x, *popt)
            pp.plot(x, yfit, '-')
            pp.title((str(r.params)))
            pp.xlabel(str(popt))

            # Add it to the dictionary
            add_result_to_dict(out_dict, speedup(fit_func, popt), 1.0, r.params)
            n_fit += 1
        except:
            pp.title(str(r.params))
            pp.xlabel("Didn't fit")

            add_result_to_dict(out_dict, 0, 0.0, r.params)

        pp.xlim(xmin=0)
        pp.ylim(ymin=0)
        pp.savefig(os.path.join(fig_out_dir, 'fit-%d.png' % i))
        log.debug("Fit result %d", i)

    return out_dict



def average_many(superrootdir='.', regexp='d-run-[0-9]+'):

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
    out_dict = average_many()
    with open('d-run.list.pickl', 'w') as f:
        pickle.dump(out_dict, f)

if __name__=="__main__":
    log.basicConfig(level=log.INFO)
    main()