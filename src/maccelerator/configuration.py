"""Configure types of runs."""

from simtk import openmm
import scipy.io


class Configuration(object):
    def __init__(self):
        self.simulator = None
        self.modeller = None
        self.adapter = None

    def get_param_grid(self):
        raise NotImplementedError


class OpenMMConfiguration(Configuration):
    def serialize(self, system, integrator, out_sys_fn, out_int_fn):
        #TODO: This code should probably go elsewhere

        with open(out_sys_fn, 'w') as f:
            f.write(openmm.XmlSerializer.serialize(system))
        with open(out_int_fn, 'w') as f:
            f.write(openmm.XmlSerializer.serialize(integrator))


class TMatConfiguration(Configuration):
    def __init__(self, tmat_fn):
        self.tmat = scipy.io.mmread(tmat_fn)
        self.tmat = self.tmat.tocsr()





