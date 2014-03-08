'''
Created on Nov 18, 2013

@author: harrigan
'''

import os
import sys

PBS_HEADER = """
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l mem=8gb
#PBS -j oe
#PBS -o {d}
#PBS -M harrigan@stanford.edu
#PBS -m ae

cd $PBS_O_WORKDIR/{d}
export OMP_NUM_THREADS=1

"""

MSM_FROM_ASSIGNMENTS_SCRIPT = """
python -u ../../../src/quantaccel/build_msm_from_assignments.py {round_i} {which} {how} &> tmatfromass-{how}-{round_i:d}.log

"""

MSM_FROM_CLUSTER_SCRIPT = """
build_msm_from_cluster.py {round_i} {how} &> tmatfromclus-{how}-{round_i:d}.log
"""

SUBMITTER_HEADER = """
PREVWD=$CWD
"""



JOBFN = "{whence}-{how}-mk6-{round_i:d}.job"
OUTFN = "{whence}-{how}-mk6-{round_i:d}.mtx"
SUBMITFN = 'submit.sh'

SUBMITTER_ENTRY = """
if [ ! -e {outfn} ]
    then
        mqsub {jobfn}
fi

"""

def do(dir_dat, qsub_header, executable, options, finder='models'):
    with open(dir_dat) as f:
        dirs = [line.strip() for line in f if len(line) > 1]

    with open (SUBMITFN, 'w') as subf:
        for d in dirs:
            nrounds = len(os.listdir(os.path.join(d, finder)))
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
        do('dirs_round_new.dat', PBS_HEADER, MSM_FROM_CLUSTER_SCRIPT, {'which':'muller', 'how':'round', 'whence':'tmatfromclus'})
    elif sys.argv[1] == 'mullerclusolt':
        do('dirs_lts.dat', PBS_HEADER, MSM_FROM_CLUSTER_SCRIPT, {'which':'muller', 'how':'percent', 'whence':'tmatfromclus'})
    elif sys.argv[1] == 'muller3':
        do('dirs.dat', PBS_HEADER, MSM_FROM_CLUSTER_SCRIPT, {'how':'rnew', 'whence': 'tmat'}, finder='trajs')
