'''
Created on Nov 18, 2013

@author: harrigan
'''

import os
import sys

PBS_HEADER = """
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -l mem=8gb
#PBS -j oe
#PBS -o {d}
#PBS -M harrigan@stanford.edu
#PBS -m bea

cd $PBS_O_WORKDIR/{d}
export OMP_NUM_THREADS=1

"""

MSM_FROM_ASSIGNMENTS_SCRIPT = """
python -u ../../../src/quantaccel/build_msm_from_assignments.py {round_i} {which} {how} &> tmatfromass-{how}-{round_i:d}.log

"""

MSM_FROM_CLUSTER_SCRIPT = """
python -u ../../../src/quantaccel/build_msm_from_cluster.py {round_i} {which} {how} &> tmatfromclus-{how}-{round_i:d}.log
"""

SUBMITTER_HEADER = """
PREVWD=$CWD
"""



JOBFN = "{whence}-{how}-mkiii-{round_i:d}.job"
OUTFN = "{whence}-{how}-mkiii-{round_i:d}.mtx"
SUBMITFN = 'submit.sh'

SUBMITTER_ENTRY = """
if [ ! -e {outfn} ]
    then
        qsub {jobfn}
        sleep 1
fi

"""

def do(dir_dat, qsub_header, executable, options):
    with open(dir_dat) as f:
        dirs = [line.strip() for line in f if len(line) > 1]

    with open (SUBMITFN, 'w') as subf:
        for d in dirs:
            nrounds = len(os.listdir(os.path.join(d, 'models')))
            print "%d rounds for %s" % (nrounds, d)

            for i in xrange(nrounds):
                jobfn = os.path.join(d, JOBFN.format(round_i=i, **options))
                outfn = os.path.join(d, OUTFN.format(round_i=i, **options))
                with open(jobfn, 'w') as f:
                    f.write(qsub_header.format(d=d))
                    f.write(executable.format(round_i=i, **options))
                subf.write(SUBMITTER_ENTRY.format(outfn=outfn, jobfn=jobfn))


if __name__ == '__main__':
    if sys.argv[1] == 'tmatass':
        do('dirs.dat', PBS_HEADER, MSM_FROM_ASSIGNMENTS_SCRIPT, {'which': 'tmat', 'how':'round', 'whence':'tmatfromass'})
    elif sys.argv[1] == 'mullerass':
        do('dirs.dat', PBS_HEADER, MSM_FROM_ASSIGNMENTS_SCRIPT, {'which': 'muller', 'how': 'round', 'whence':'tmatfromass'})
    elif sys.argv[1] == 'tmatolt':
        do('oltdirs.dat', PBS_HEADER, MSM_FROM_ASSIGNMENTS_SCRIPT, {'which': 'tmat', 'how': 'percent', 'whence':'tmatfromass'})
    elif sys.argv[1] == 'mullerclus':
        do('dirs.dat', PBS_HEADER, MSM_FROM_CLUSTER_SCRIPT, {'which':'muller', 'how':'round', 'whence':'tmatfromclus'})
    elif sys.argv[1] == 'mullerclusolt':
        do('oltdirs.dat', PBS_HEADER, MSM_FROM_CLUSTER_SCRIPT, {'which':'muller', 'how':'percent', 'whence':'tmatfromclus'})
