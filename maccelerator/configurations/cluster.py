__author__ = 'harrigan'

from .base import ClusterConfig

_PBS_SCRIPT = """#!/bin/bash
#PBS -l walltime={hours}:00:00
#PBS -j {join_output}
#PBS -l nodes=1:ppn={ppn}

cd $PBS_O_WORKDIR
{run_line}
"""

_SLURM_SCRIPT = """#!/bin/bash
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --job-name=accelerate
#SBATCH --output=accelerate-%j.sout
#SBATCH --time={days}-{hours}:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={ppn}

cd $SLURM_SUBMIT_DIR
{run_line}
"""

_SERIAL_RUN = """
for i in {{1..{n_copy}}}
do
    python {script_fn} -i $i
done
"""

_GNU_PARALLEL_RUN = """
seq {n_copy} | parallel -j {ppn} python {script_fn} -i {{}}"""


class PBSCluster(ClusterConfig):
    def __init__(self, n_copy, parallel='parallel', ppn=16, hours=72,
                 join_output='oe'):
        self.n_copy = n_copy
        self.ppn = ppn
        self.hours = hours
        self.join_output = join_output
        self.parallel = str(parallel).lower()

    @property
    def job_script_ext(self):
        return 'job'

    def make_job_script(self, py_fn):

        # Choose parallel type
        if self.parallel == 'parallel':
            fmt_dict_r = dict(n_copy=self.n_copy, ppn=self.ppn, script_fn=py_fn)
            run_line = _GNU_PARALLEL_RUN.format(**fmt_dict_r)
        elif self.parallel == 'serial':
            fmt_dict_r = dict(n_copy=self.n_copy, script_fn=py_fn)
            run_line = _SERIAL_RUN.format(**fmt_dict_r)
        else:
            raise ValueError('Invalid parallel type: %s', self.parallel)

        # Return the job script
        fmt_dict_j = dict(hours=self.hours, join_output=self.join_output,
                          ppn=self.ppn, run_line=run_line)
        return _PBS_SCRIPT.format(**fmt_dict_j)


class SlurmCluster(ClusterConfig):
    def __init__(self, n_copy, parallel='parallel', ppn=16, days=2, hours=0):
        self.n_copy = n_copy
        self.ppn = ppn
        self.days = days
        self.hours = hours
        self.parallel = str(parallel).lower()

    @property
    def job_script_ext(self):
        return 'sbatch'

    def make_job_script(self, py_fn):
        # Choose parallel type
        if self.parallel == 'parallel':
            fmt_dict_r = dict(n_copy=self.n_copy, ppn=self.ppn, script_fn=py_fn)
            run_line = _GNU_PARALLEL_RUN.format(**fmt_dict_r)
        elif self.parallel == 'serial':
            fmt_dict_r = dict(n_copy=self.n_copy, script_fn=py_fn)
            run_line = _SERIAL_RUN.format(**fmt_dict_r)
        else:
            raise ValueError('Invalid parallel type: %s', self.parallel)

        # Return job script
        job_fmt_dict = dict(hours=self.hours, days=self.days, ppn=self.ppn,
                            run_line=run_line)
        return _SLURM_SCRIPT.format(**job_fmt_dict)

