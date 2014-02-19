# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:19:46 2014

@author: harrigan
"""
import matplotlib
matplotlib.use('Agg')

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
        ma1 = re.match(r"^result-h-runcopy-[0-9]+-[0-9]+\.pickl$", fn)
        ma2 = re.match(r"^result-ho-[0-9]+-[0-9]+\.pickl$", fn)
        if ma1 or ma2:
            with open(os.path.join(rootdir, fn)) as f:
                r = pickle.load(f)
                r.params['adaptive'] = True if ma1 else False
                results.append(fix_result(r))
    return results


def fix_result(r):
    """Process the results to make them fittable

    e.g. make real, remove outliers.
    """
    if not np.all(np.isreal(r.errors)):
        raise Exception("Unreal!")

    #r.errors = np.delete(r.errors, np.where(r.errors[:, 1] > 1.0)[0], 0)
    return r


def cliff_find(errors, highc, lowc, percent=False):
    """Find the x val that is the cliff."""

#    inds = np.lexsort((errors[:,1], errors[:,0]))
#    errors = errors[inds]

    vals = errors[:, 1]
    highscore = 0.0
    higharg = -1
    for i in range(len(vals)+1):
        left = vals[:i]
        n_left = len(left)
        right = vals[i:]
        n_right = len(right)

        # Should we sum percentages or totals?
        n_high = np.sum(left > highc)
        n_low = np.sum(right < lowc)

        if percent:
            # Percentage
            if n_high == 0 and n_left == 0:
                leftscore = 1.0
            else:
                leftscore = n_high / n_left

            if n_low == 0 and n_right == 0:
                rightscore = 1.0
            else: rightscore = n_low / n_right

            score = 0.5 * (leftscore + rightscore)
        else:
            # Total
            score = n_high + n_low

        if score > highscore:
            highscore = score
            higharg = i


    if higharg < len(vals) and higharg > 0:
        cliff_x = (errors[higharg, 0] + errors[higharg - 1, 0]) / 2.0
    elif higharg == len(vals):
        cliff_x = np.inf
        # Bias = bias of the last element, which is high
        bias = errors[-1,1]
        return cliff_x, [0, 0, bias], np.linspace(0,np.max(errors[:,0]), 2), np.array([bias]*2)
    elif higharg == 0:
        cliff_x = errors[0, 0]

    # Try to fit
    popt, pcov = scipy.optimize.curve_fit(
        offset_exp, errors[higharg:, 0], errors[higharg:, 1], p0=[100,1,0])

    x = np.linspace(cliff_x, np.max(errors[:,0]), 100)
    y = offset_exp(x, *popt)

    tau, c, off = popt
    return cliff_x, popt, x, y

def ifunc(x, a):
    """Inverse."""
    return a / x


def ifunc_off(x, a, offset):
    """Inverse with a y-offset."""
    return a / x - offset


def simple_exp(x, tau):
    """Exponential"""
    return np.exp(-x / tau)

def offset_exp(x, tau, c, off):
    return c*np.exp(-x / tau) + off

def bendy_exp(x, tau, n):
    """Exponential where the power of the exponent can vary."""
    return np.exp(- (x / tau) ** n)

def offset_bendy_exp(x, tau, n, off):
    return (1-off) * np.exp(- (x / tau) ** n) + off

def expifunc(x, a, n):
    """Inverse with tunable exponent."""
    return a / (x ** n)

def double_exp(x, c, d, tau1, tau2):
    return c * np.exp(-x/tau1) + d * np.exp(-x/tau2)

def sigmoid(t, tau, t0, c, off):
    return ((c-off)/(1+np.exp((t-t0)/tau))) + off

def two_exps(t, tau1, tau2, off1, off2, a1, a2, t0):

    right = (a2 - off2) * np.exp(-t/tau2) + off2
    left = (a1 - off1) * np.exp(-t/tau1) + off1
    mright = 0.5 * (np.sign(t-t0)+1)
    mleft = 0.5 * (np.sign(-t+t0) + 1)
    return mleft * left + mright * right

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
    elif fit_func == offset_bendy_exp:
        # TODO: (maybe) put something better here
        tau, n, off = popt
        time = tau * (-np.log(err_bits))**(1.0/n)
        return scale / tau
    elif fit_func == sigmoid:
        tau, t0, c, off = popt
        return scale / t0
    elif fit_func == two_exps:
        tau1, tau2, off1, off2, a1, a2, t0 = popt
        return scale / t0
    elif fit_func == cliff_find:
        cliff_x, tau, c, off = popt
        return scale / cliff_x
    else:
        raise Exception("No speedup calc for this fit func.")


def _params_to_str(params, get_rid_of=['lag_time', 'runcopy']):
    return ', '.join(["%s: %s" % (key, val)
                        for (key, val) in params.items()
                        if key not in get_rid_of])

def fit_results(results, fig_out_dir, fit_func=cliff_find, p0=None, show=False):
    """Fit all the results from a particular runcopy (directory of results)
    """
    n_fit = 0

    if p0 is None:
        if fit_func == expifunc:
            p0 = [1000., 1.]
        elif fit_func == simple_exp:
            p0 = [1000.]
        elif fit_func == bendy_exp:
            p0 = [1000., 1.]
        elif fit_func == offset_bendy_exp:
            p0 = [100., 0.5, 0.8]
        elif fit_func == sigmoid:
            p0 = [1, 500, 1, 0.1]
        elif fit_func == two_exps:
            p0 = [1000, 1000, 0.8, 0.1, 1.0, 0.2, 1000]
        elif fit_func == cliff_find:
            p0 = [0.6, 0.4]
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

            if fit_func == cliff_find:
                highc, lowc = p0
                cliff_x, popt, x, y = cliff_find(r.errors, highc, lowc)
                pp.hlines(p0, 0.0, np.max(r.errors[:,0]))
                pp.vlines(cliff_x, 0, 1.0)

                pp.plot(x, y, 'r-')
                tau, c, off = popt
                quality = 1.0 - off
                popt = np.append([cliff_x], popt)

            else:
                # Try to fit
                popt, pcov = scipy.optimize.curve_fit(
                    fit_func, r.errors[:, 0], r.errors[:, 1], p0=p0)
                # Make fit curve
                if fit_func in [expifunc]:
                    # For those that blow up at 0
                    x = np.linspace(np.min(r.errors[:, 0]), np.max(r.errors[:, 0]))
                else:
                    x = np.linspace(0, np.max(r.errors[:,0]), num=200)
                yfit = fit_func(x, *popt)
                pp.plot(x, yfit, '-')

                # If the parameters do not change, it probably didn't fit
                if np.all(np.abs(popt - p0) < 1e-6):
                    raise Exception("Same as p0")

                # We don't want any negative params either
                if np.any(popt<0):
                    raise Exception("Negative params")

                quality = 1.0

            speedup = calc_speedup(fit_func, popt)

            # Labels
            pp.title(_params_to_str(r.params))
            pp.xlabel(', '.join(["%f" % po for po in popt]))
            pp.suptitle("Speedup: %d, Quality: %f" % (int(speedup), quality))

            # Add it to the dictionary
            add_result_to_dict(
                out_dict, speedup,
                quality, r.params)
            n_fit += 1
        except Exception as e:
            pp.title(_params_to_str(r.params))
            pp.suptitle("Didn't fit: %s" % e.args)

            add_result_to_dict(out_dict, 0, 0.0, r.params)

        pp.xlim(xmin=0)
        ymin, ymax = pp.ylim()
        pp.ylim(ymin=min(ymin, 0.0), ymax=max(ymax,1.0))
        pp.savefig(os.path.join(fig_out_dir, 'fit-%d.png' % i))
        if show: pp.show()
        log.debug("Fit result %d", i)

    log.info("Found fits for %d out of %i", n_fit, len(results))

    return out_dict


def average_many(superrootdir='.', regexp='^h-run-[0-9]+$', test=False):
    """Fit over a set of runcopys (directories)."""

    if test:
        regexp='^h-run-1$'

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
    with open('h-run.list.pickl', 'w') as f:
        pickle.dump(out_dict, f)

if __name__ == "__main__":
    log.basicConfig(level=log.INFO)
    main()
