"""Keep track of the directory structure and filenames.

For adaptive runs, we have to deal with a relatively elaborate directory
structure. Hardcoded configuration-independent path information should be
contained here.

If there are filenames which may differ depending on the configuration
(e.g. traj.dcd vs. traj.h5), they should be properties of the relevant
object (e.g. traj_fn is a property of Simulator)
"""

__author__ = 'harrigan'

from os.path import join as pjoin
import os
import logging as log


class FileStructure():
    """Class to deal with folder structure and filenames.

    ./
        trajs/
            round-0/
                traj-0.file
                traj-1.file
            round-1/
                ...
            ...
        sstates/
            sstate-0.file
            sstate-1.file
            ...
        msms/
            msm.adapt-0.file
            msm.conv-0.file
            msm.adapt-1.file
            msm.conv-1.file
            ...
        figs/
            plot-0000.png
            plot-0001.png
            ...

    Note: In the future, if there is desire to create different folder
    structures, this class can be subclassed similar to the other objects.
    """

    # Base names
    traj_base = 'trajs'
    sstate_base = 'sstates'
    msms_base = 'msms'
    figs_base = 'figs'

    trajround_fmt = 'round-{round_i}'

    def __init__(self):
        self.traj_dir = ''
        self.sstate_dir = ''
        self.msms_dir = ''
        self.figs_dir = ''
        self.rundir = ''

    def make_directories(self, rundir):
        """Set up directories for a run."""
        self.rundir = rundir
        self.traj_dir = pjoin(rundir, self.traj_base)
        self.sstate_dir = pjoin(rundir, self.sstate_base)
        self.msms_dir = pjoin(rundir, self.msms_base)
        self.figs_dir = pjoin(rundir, self.figs_base)
        try:
            os.mkdir(self.traj_dir)
            os.mkdir(self.sstate_dir)
            os.mkdir(self.msms_dir)
            os.mkdir(self.figs_dir)
        except OSError as e:
            log.warning('Skipping %s (%s)', rundir, e)
            return False
        return True

    def make_trajround(self, round_i):
        """Make directory for trajectories for a round"""
        trajround_dir = self.trajround_fmt.format(round_i=round_i)
        trajround_dir = pjoin(self.traj_dir, trajround_dir)
        os.mkdir(trajround_dir)
        return trajround_dir
