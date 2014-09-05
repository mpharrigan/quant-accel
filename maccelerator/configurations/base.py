"""Configure types of runs."""

from ..files import FileStructure
import pickle
import mdtraj
import re
import inspect

_GENERAL_TEMPLATE = """from maccelerator import *
import itertools
import argparse
import logging

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=int, default=0)
args = parser.parse_args()
logging.basicConfig(level=logging.INFO, filename='copy-{{}}.log'.format(args.i))

config = {config_class}()
{other_config}
{defaults}

config.apply_configuration()
with {grid_manager}(config, run_id=args.i) as grid:
    grid.grid()
"""


class ClusterConfig:
    """Base class for configuring supercomputer clusters."""

    @property
    def job_script_ext(self):
        raise NotImplementedError

    def make_job_script(self, py_fn):
        raise NotImplementedError


class Configuration(object):
    def __init__(self):
        self.simulator = None
        self.modeller = None
        self.convchecker = None
        self.adapter = None

        self.simulator_class = None
        self.modeller_class = None
        self.convchecker_class = None
        self.adapter_class = None

        # Right now, all configurations use the same file structure.
        # If we want to introduce new file structure in the future,
        # move this object out of the base class and treat it like all
        # the other objects.
        self.file = FileStructure(self)

        self.defaults()

    def apply_configuration(self):
        # TODO: Apply here after ducking every object to take "config"
        raise NotImplementedError

    def get_param_grid(self, run_id):
        raise NotImplementedError

    @classmethod
    def get_template(cls, grid_manager_name):
        """Use the method self.defaults to write a configuration."""

        # Get code from defaults method
        code, _ = inspect.getsourcelines(cls.defaults)

        # Remove two levels of indentation and comment everything
        code_parsed = [re.sub(r'^        ', '# ', c) for c in code[1:]]

        # Parse and format
        code_parsed = ''.join(code_parsed)
        return _GENERAL_TEMPLATE.format(grid_manager=grid_manager_name,
                                        defaults=code_parsed,
                                        config_class=cls.__name__,
                                        other_config=cls._other_config())

    @classmethod
    def _other_config(cls):
        """Overwrite this function to add more things to the config script."""
        return ""

    @property
    def seed_state(self):
        return self.adapter.seed_states

    @property
    def simulate(self):
        return self.simulator.simulate

    @property
    def model(self):
        return self.modeller.model

    @property
    def adapt(self):
        return self.adapter.adapt

    @property
    def check_convergence(self):
        return self.convchecker.check_convergence

    def defaults(config):
        raise NotImplementedError

    def __getstate__(self):
        state = dict(self.__dict__)
        try:
            del state['get_param_grid']
        except KeyError:
            pass
        return state


class OpenMMConfiguration(Configuration):
    pass


class TMatConfiguration(Configuration):
    def __init__(self, ref_msm_fn, centers_fn):
        super().__init__()

        # Load reference MSM
        with open(ref_msm_fn, 'rb') as ref_msm_f:
            self.ref_msm = pickle.load(ref_msm_f)

        # Load cluster centers for visualization
        self.centers = mdtraj.io.loadh(centers_fn, 'cluster_centers')




