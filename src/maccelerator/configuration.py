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
    pass




