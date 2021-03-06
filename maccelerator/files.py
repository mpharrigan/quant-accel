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
import logging
import pickle

from pkg_resources import resource_filename


log = logging.getLogger(__name__)


def get_fn(fn):
    return resource_filename('maccelerator', 'reference/{}'.format(fn))


class ConfigUnpickler(pickle.Unpickler):
    """Deal with custom subclasses."""

    def find_class(self, module, name):
        """If the main script defines a subclass whose name begins with My:
        load the class whose name lacks the 'My'
        """
        if module == '__main__' and 'My' in name and name[:2] == 'My':
            return super().find_class('maccelerator', name[2:])
        else:
            return super().find_class(module, name)


def special_pickle_load(fn):
    """Wrapper to use ConfigUnpickler

    :param fn: Filename. We take care of opening the file.
    """
    with open(fn, 'rb') as f:
        return ConfigUnpickler(f).load()


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
    convs_base = 'convs'

    trajround_fmt = 'round-{round_i}'

    def __init__(self, config):
        self.traj_dir = ''
        self.sstate_dir = ''
        self.msms_dir = ''
        self.figs_dir = ''
        self.convs_dir = ''
        self.rundir = ''
        self.config = config

    def make_directories(self, rundir):
        """Set up directories for a run."""
        self.rundir = rundir
        self.traj_dir = pjoin(rundir, self.traj_base)
        self.sstate_dir = pjoin(rundir, self.sstate_base)
        self.msms_dir = pjoin(rundir, self.msms_base)
        self.figs_dir = pjoin(rundir, self.figs_base)
        self.convs_dir = pjoin(rundir, self.convs_base)
        try:
            os.mkdir(self.traj_dir)
            os.mkdir(self.sstate_dir)
            os.mkdir(self.msms_dir)
            os.mkdir(self.figs_dir)
            os.mkdir(self.convs_dir)
        except OSError as e:
            log.warning('Skipping %s (%s)', rundir, e)
            return False
        return True

    def make_traj_fns(self, round_i, traj_is):
        """Make directory for trajectories for a round

        :param round_i: Round index for directory name
        :param traj_is: Iterable over trajectory indices

        """
        trajround_dir = self.trajround_fmt.format(round_i=round_i)
        trajround_dir = pjoin(self.traj_dir, trajround_dir)
        os.mkdir(trajround_dir)
        return [
            pjoin(trajround_dir, self.config.simulator.trajfn.format(traj_i=i))
            for i in traj_is]

    def sstate_fn(self, round_i):
        """Return a starting-state filename."""
        return pjoin(self.config.file.sstate_dir,
                     self.config.adapter.sstatefn.format(round_i=round_i + 1))

    def plot_fn(self, round_i, rel=False):
        """Return a plot filename."""
        relname = self.config.convchecker.plotfn.format(round_i=round_i)
        if rel:
            return pjoin(self.figs_base, relname)
        else:
            return pjoin(self.figs_dir, relname)

    def model_fn(self, round_i):
        """Return a model filename for each round."""
        return pjoin(self.msms_dir,
                     self.config.modeller.modelfn.format(round_i=round_i))

    def conv_fn(self, round_i, rel=False):
        """Return a convergence filename for each round."""
        relname = self.config.convchecker.convfn.format(round_i=round_i)
        if rel:
            return pjoin(self.convs_base, relname)
        else:
            return pjoin(self.convs_dir, relname)

    def param_fn(self, param):
        """Return a filename to save param info."""
        return pjoin(self.rundir, param.paramfn.format())

    def run_fn(self, run):
        """Return a filename to save a run object."""
        return pjoin(self.rundir, run.runfn.format())

