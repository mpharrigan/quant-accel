"""Configure types of runs."""


class Configuration(object):
    def __init__(self):
        self.simulator = None
        self.modeler = None

    @staticmethod
    def get_configs():
        raise NotImplementedError


class OpenMMConfiguration(Configuration):
    pass


class TMatConfiguration(Configuration):
    pass





