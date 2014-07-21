"""Deal with parameters."""

from os.path import join as pjoin
import os
import logging

log = logging.getLogger()


class AdaptiveParams(object):
    TRAJ = 'trajs'
    SSTATE = 'sstates'
    MSMS = 'msms'
    FIGS = 'figs'

    TRAJROUND = 'round-{round_i}'


    def __init__(self, spt, tpr, run_id=0):
        self.tpr = tpr
        self.spt = spt
        self.run_id = run_id
        self.rundir = None
        self.traj_dir = ""
        self.sstate_dir = ""
        self.msms_dir = ""
        self.figs_dir = ""


    @property
    def post_converge(self):
        raise NotImplementedError

    @property
    def adapt_lt(self):
        raise NotImplementedError

    @property
    def build_lt(self):
        raise NotImplementedError

    @property
    def dirname(self):
        return "blt-{build_lt}_alt-{adapt_lt}_spt-{spt}_tpr-{tpr}".format(
            build_lt=self.build_lt, adapt_lt=self.adapt_lt, spt=self.spt,
            tpr=self.tpr
        )

    @property
    def pretty_desc(self):
        """A description e.g. for plot titles."""
        return self.dirname

    # TODO: Should all this directory stuff really be here? Params should be
    # TODO: things that change over the grid. This stuff shouldn't change!
    # TODO: Maybe put in Configuration or make a new FileStructure class
    def make_directories(self):
        """Set up directories for a run."""
        rundir = self.rundir
        self.traj_dir = pjoin(rundir, self.TRAJ)
        self.sstate_dir = pjoin(rundir, self.SSTATE)
        self.msms_dir = pjoin(rundir, self.MSMS)
        self.figs_dir = pjoin(rundir, self.FIGS)
        try:
            os.mkdir(self.traj_dir)
            os.mkdir(self.sstate_dir)
            os.mkdir(self.msms_dir)
            os.mkdir(self.figs_dir)
        except OSError as e:
            log.warning('Skipping %s (%s)', rundir, e)
            return '', False
        return self.traj_dir, True

    def make_trajround(self, round_i):
        """Make directory for trajectories for a round"""
        trajround_dir = pjoin(self.traj_dir,
                              self.TRAJROUND.format(round_i=round_i))
        trajround_dir = os.path.abspath(trajround_dir)
        os.mkdir(trajround_dir)
        return trajround_dir
