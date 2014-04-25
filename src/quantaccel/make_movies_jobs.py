"""Make jobs for making movies."""

__author__ = 'harrigan'

import logging as log

log.basicConfig(level=log.INFO)

PBS_HEADER = """
#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l mem=8gb
#PBS -j oe
#PBS -M harrigan@stanford.edu
#PBS -m ae

: ${{PBS_O_WORKDIR="$PWD"}}
cd $PBS_O_WORKDIR
export PBS_O_WORKDIR
export OMP_NUM_THREADS=1
"""

MOVIE_JOB = PBS_HEADER + """
movietype="{movie_type}"

while read p
do
    pdir=`dirname $p`
    echo $pdir
    cd $pdir
    python -m quantaccel.centroid_movie ./msms 7 -how rnew --movietype $movietype
    python -m quantaccel.centroid_movie ./msms 7 -how pnew --movietype $movietype
    cd $PBS_O_WORKDIR
done < param_locs-{job_i}.dat

echo "Finished"
"""

SUB_LINE = """mqsub {job_fn}"""

import collections
import os
import stat
import argparse


def make_jobs(param_locs_fn='../param_locs.dat',
              param_locs_out_fn='param_locs-{job_i}.dat',
              job_fn='make-movies-{movie_type}-{job_i}.job',
              movies_per_job=50,
              movie_type='projection-pop'):
    """Take a long list of locations and divvy it up.

    :param param_locs_fn: Where to find the complete list of directories
    :param param_locs_out_fn: Where to put our divvy'd files
    :param job_fn: The format string for job filenames
    :param movies_per_job: How many movies to do in one job
    :param movie_type: What type of movie to make. This will do both pnew and
                        rnew
    """

    # Data structure to collect out files
    divvy_dat = collections.defaultdict(list)

    # Assign to data structure
    with open(param_locs_fn) as param_locs_f:
        for i, line in enumerate(param_locs_f):
            job_i = i // movies_per_job
            divvy_dat[job_i].append(line)

    # Write files
    submit_lines = []
    for i, dat_lines in divvy_dat.items():
        cur_job_fn = job_fn.format(job_i=i, movie_type=movie_type)
        with open(cur_job_fn, 'w') as job_f:
            job_f.write(MOVIE_JOB.format(movie_type=movie_type, job_i=i))

        with open(param_locs_out_fn.format(job_i=i), 'w') as dat_f:
            dat_f.writelines(dat_lines)

        submit_lines += [SUB_LINE.format(job_fn=cur_job_fn)]

    with open('submit.sh', 'w') as sub_f:
        sub_f.write('\n'.join(submit_lines))

    # Make executable
    st = os.stat(os.path.join('.', 'submit.sh'))
    os.chmod(os.path.join('.', 'submit.sh'), st.st_mode | stat.S_IEXEC)


def parse():
    """Use argparse for command line parsing."""
    parser = argparse.ArgumentParser(description="Make jobs for making movies",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('movie_type',
                        help="""What type of movie to make: [projection-pop,
                        centroid-it, num-states]""")

    args = parser.parse_args()
    make_jobs(movie_type=args.movie_type)


if __name__ == "__main__":
    make_jobs()


