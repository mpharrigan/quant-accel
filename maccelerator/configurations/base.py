"""Configure types of runs."""

from ..files import FileStructure

_GENERAL_TEMPLATE = """import maccelerator as maccel
import itertools
import argparse
import logging

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=int, default=0)
args = parser.parse_args()
logging.basicConfig(level=logging.INFO, filename='copy-{{}}.log'.format(args.i))

# Important: any subclasses defined here must be named `My{{name of superclass}}`
{other_config}

with maccel.{grid_manager}(MyAlanineConfiguration(), run_id=args.i) as grid:
    grid.grid()
"""


class ClusterConfig:
    pass


class Configuration(object):
    def __init__(self):
        self.simulator = None
        self.modeller = None
        self.convchecker = None
        self.adapter = None

        # Right now, all configurations use the same file structure.
        # If we want to introduce new file structure in the future,
        # move this object out of the base class and treat it like all
        # the other objects.
        self.file = FileStructure(self)

    def get_param_grid(self, run_id):
        raise NotImplementedError

    @classmethod
    def get_template(cls, grid_manager_name):
        return _GENERAL_TEMPLATE

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


class OpenMMConfiguration(Configuration):
    pass


class TMatConfiguration(Configuration):
    pass




