"""Configure types of runs."""

from simtk import openmm
import scipy.io


class Configuration(object):
    def __init__(self):
        self.simulator = None
        self.modeller = None
        self.convchecker = None
        self.adapter = None

    def get_param_grid(self):
        raise NotImplementedError


class OpenMMConfiguration(Configuration):
    pass


class TMatConfiguration(Configuration):
    def __init__(self, tmat_fn):
        super().__init__()
        self.tmat = scipy.io.mmread(tmat_fn)
        self.tmat = self.tmat.tocsr()





