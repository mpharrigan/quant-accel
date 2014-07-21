"""Configure types of runs."""

from ..files import FileStructure

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
        self.file = FileStructure()

    def get_param_grid(self):
        raise NotImplementedError


class OpenMMConfiguration(Configuration):
    pass


class TMatConfiguration(Configuration):
    pass




