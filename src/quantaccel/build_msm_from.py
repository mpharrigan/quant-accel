"""This function calls either build_msm_from_assignment or
build_msm_from_cluster."""

__author__ = 'harrigan'

import sys
from . import build_msm_from_assignments, build_msm_from_cluster

if __name__ == "__main__":
    _, system_type, round_i, how = sys.argv
    round_i = int(round_i)
    if system_type == 'muller':
        build_msm_from_cluster.do(round_i, how)
    elif system_type == 'tmat':
        build_msm_from_assignments.do(round_i, how)
    else:
        raise ValueError('Invalid system type')
