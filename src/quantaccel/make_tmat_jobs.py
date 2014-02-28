import sys
import numpy as np
import os

JOB_SCRIPT = """

#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l mem=8gb
#PBS -j oe
#PBS -o .
#PBS -M harrigan@stanford.edu
#PBS -m a

# We need 80

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1

RUNCOPY={runcopy}

python ../../src/quantaccel/tmat_simulation.py 0 {job_i} $RUNCOPY &> tmat-{job_i}.log
"""

STARTSTATE_FN = 'starting_state.int'
UNFOLD_FN = '../unfolded.dat'

def main(runcopy, num_permutes=80):

    # Pick a starting state for all of them
    if not os.path.exists(STARTSTATE_FN):
        unfolded = np.loadtxt(UNFOLD_FN, dtype=int)

        # start_state = np.random.randint(0, 10000)
        start_state = np.random.choice(unfolded)
        with open(STARTSTATE_FN, 'w') as f:
            f.write("{}".format(start_state))

    # Write all the job files
    for job_i in range(num_permutes):
        with open('tmat-{job_i}.job'.format(job_i=job_i), 'w') as f:
            f.write(JOB_SCRIPT.format(job_i=job_i, runcopy=runcopy))


if __name__ == "__main__":
    main(runcopy=sys.argv[1])
