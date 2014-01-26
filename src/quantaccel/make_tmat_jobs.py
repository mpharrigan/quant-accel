import sys

JOB_SCRIPT = """

#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -l mem=8gb
#PBS -j oe
#PBS -o .
#PBS -M harrigan@stanford.edu
#PBS -m bea

# We need 60

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1

RUNCOPY={runcopy}

python ../../src/quantaccel/tmat_simulation.py 0 {job_i} $RUNCOPY &> tmat-{job_i}.log
"""

def main(runcopy, num_permutes=60):
    for job_i in range(num_permutes):
        with open('tmat-{job_i}.job'.format(job_i=job_i), 'w') as f:
            f.write(JOB_SCRIPT.format(job_i=job_i, runcopy=runcopy))


if __name__ == "__main__":
    main(runcopy=sys.argv[1])