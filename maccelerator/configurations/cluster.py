__author__ = 'harrigan'

from .base import ClusterConfig

_PBS_SCRIPT = """#!/bin/bash
#PBS -l walltime={hours}:00:00
#PBS -j {join_output}
#PBS -l nodes=1:ppn={ppn}

cd $PBS_O_WORKDIR
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


