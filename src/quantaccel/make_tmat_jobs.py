
import numpy as np
import os
import argparse
import re

JOB_SCRIPT = """

#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -l mem=8gb
#PBS -j oe
#PBS -o .
#PBS -M harrigan@stanford.edu
#PBS -m a


cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1


python -m quantaccel.tmat_simulation {job_i} --runcopy {runcopy} --adapt {adapt} &> tmat-{job_i}.log
"""

QSUB = """
if [ ! -e {outfn} ]
    then
        mqsub {jobfn}
fi
"""


STARTSTATE_FN = 'starting_state.int'
UNFOLD_FN = '../unfolded.dat'
NADAPTIVE = 30
NCONTROL = 5

def write_jobs(runcopy, adapt):
    """Write jobs

    Number given by global variables.
    """

    # Pick a starting state for all of them
    if not os.path.exists(STARTSTATE_FN):
        unfolded = np.loadtxt(UNFOLD_FN, dtype=int)

        # start_state = np.random.randint(0, 10000)
        start_state = np.random.choice(unfolded)
        with open(STARTSTATE_FN, 'w') as f:
            f.write("{}".format(start_state))

    # Write all the job files
    for job_i in range(NADAPTIVE):
        with open('tmat-{adapt}-{job_i}.job'.format(job_i=job_i, adapt=adapt), 'w') as f:
            f.write(JOB_SCRIPT.format(job_i=job_i, runcopy=runcopy, adapt=adapt))

    for job_i in range(NCONTROL):
        with open('tmat-non-{job_i}.job'.format(job_i=job_i), 'w') as f:
            f.write(JOB_SCRIPT.format(job_i=job_i, runcopy=runcopy, adapt='non'))



def submit_jobs(runcopy):
    """Submit jobs found via regex."""
    dirlist = os.listdir('.')

    with open('submitter.sh', 'w') as f:
        for fn in dirlist:
            match = re.match(r'^tmat-(.+)-([0-9]+)\.job$', fn)
            if match:
                adapt = match.group(1)
                i = match.group(2)                
                outfn = 'result-{adapt}-{runcopy}-{i}.pickl'.format(
                                adapt=adapt, runcopy=runcopy, i=i)
                f.write(QSUB.format(jobfn=fn, outfn=outfn))

def parse():
    parser = argparse.ArgumentParser(description="Make tmat jobs",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('runcopy', type=int,
                        help='ID of this copy of run')
    parser.add_argument('--adapt',
                        help="What means of doing adaptive run")
    args = parser.parse_args()


    write_jobs(args.runcopy, args.adapt)
    submit_jobs(args.runcopy)


if __name__ == "__main__":
    parse()
