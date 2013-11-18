'''
Created on Nov 18, 2013

@author: harrigan
'''

PBS_HEADER = """
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -l mem=8gb
#PBS -j oe
#PBS -o $PBS_O_WORKDIR
#PBS -M harrigan@stanford.edu
#PBS -m bea
#PBS -N {jobname}

cd $PBS_O_WORKDIR

"""

MSM_FROM_ASSIGNMENTS_SCRIPT = """
python -u ../../../src/quantaccel/build_msm_from_assignments.py {round_i} {which} &> tmatfromass-round-{round_i:d}.log
"""

def do(dir_dat):
    with open(dir_dat) as f:
        dirs = [line.strip() for line in f if len(line)>1]
        

if __name__ == '__main__':
    pass