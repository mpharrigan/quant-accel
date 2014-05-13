__author__ = 'harrigan'

from .run import MAccelRun


class MAccelGrid(object):
    """Run many parameter configurations."""

    def __init__(self, configuration):
        self.config = configuration

    def grid(self):
        for params in self.config.get_param_grid():
            run = MAccelRun(self.config, params)
            run.run()
