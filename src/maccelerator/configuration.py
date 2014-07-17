"""Configure types of runs."""

from simtk import openmm
import scipy.io
import pickle


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
    def __init__(self, ref_msm):
        super().__init__()

        self.ref_msm = ref_msm
        self.tmat = self.ref_msm.transmat_.tocsr()




