'''
Created on Apr 7, 2014

@author: harrigan
'''

import argparse
import logging
import os
import stat

from quantaccel import build_msm_from_cluster, build_msm_from_assignments


log = logging.getLogger()
log.setLevel(logging.INFO)

PBS_HEADER = """#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime={hours}:00:00
#PBS -l mem=8gb
#PBS -j oe
#PBS -M harrigan@stanford.edu
#PBS -m ae

: ${{PBS_O_WORKDIR="$PWD"}}
cd $PBS_O_WORKDIR
export PBS_O_WORKDIR
export OMP_NUM_THREADS=1

"""

JOB = PBS_HEADER + """

python -m quantaccel.subdivide_one_round -st {system_type} subdivide --dirname {abspath}
"""


def subdivide(dirname, system_type):
    """Build an msm at points for a run at dirname.
    """
    os.chdir(dirname)
    n_points = build_msm_from_cluster.NPOINTS

    if system_type == 'muller':
        buildmsm_func = build_msm_from_cluster.do
    elif system_type == 'tmat':
        buildmsm_func = build_msm_from_assignments.do
    else:
        raise ValueError('Invalid system type')

    for percent_i in range(n_points):
        log.info("Doing point %d/%d", percent_i, n_points)
        buildmsm_func(percent_i, 'pnew')


def subdivide_func(args):
    """Entry point for argparse."""
    return subdivide(args.dirname, args.system_type)


def set_up_jobs(base_dir, dirlist_fn, system_type, submitter_fn='submit.sh'):
    """Write one job per configuration

    Each job will do all NPOINTS points
    """
    submitter_lines = []

    with open(dirlist_fn) as f:
        for line in f:
            line = line.strip()
            abspath = os.path.abspath(os.path.join(base_dir, line))

            # Make job filename
            job_fn = line[line.find('runcopy'):]
            job_fn = job_fn.replace('.', '').replace('/', '_')

            with open("{}.job".format(job_fn), 'w') as job_f:
                job_f.write(JOB.format(abspath=abspath,
                                       hours=6,
                                       system_type=system_type))

            submitter_lines += ['mqsub {}.job'.format(job_fn)]

    with open(submitter_fn, 'w') as sub_f:
        sub_f.write('\n'.join(submitter_lines))

    # Make executable
    st = os.stat(submitter_fn)
    os.chmod(submitter_fn, st.st_mode | stat.S_IEXEC)


def set_up_jobs_func(args):
    """Entry point for argparse."""
    set_up_jobs(args.basedir, args.dirlist_fn, args.system_type)


def parse():
    """Parse arguments."""
    parser = argparse.ArgumentParser(
        description="Take simulations that converged in one round and subdivide them",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--basedir',
                        help="Directory is relative to this",
                        default="./")

    parser.add_argument('--system_type', '-st',
                        help="Type of system: muller, tmat", default='muller')

    sp = parser.add_subparsers()
    setup_p = sp.add_parser('setup',
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    setup_p.set_defaults(func=set_up_jobs_func)
    setup_p.add_argument('--dirlist_fn',
                         help="File that lists all directories to do",
                         default="converged-in-one.dat")

    subdivide_p = sp.add_parser('subdivide',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subdivide_p.add_argument('--dirname',
                             help="Relative directory", required=True)

    subdivide_p.set_defaults(func=subdivide_func)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    parse()
