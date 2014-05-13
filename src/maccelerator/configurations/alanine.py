from ..simulate import TMatSimulator
from ..model import TMatModeller
from ..configuration import TMatConfiguration
from ..param import AdaptiveParams


class AlanineSimulator(TMatSimulator):
    def __init__(self):
        # TODO: Pass in arguments
        super().__init__()

class AlanineModeller(TMatModeller):
    pass


class AlanineParams(AdaptiveParams):
    pass


class AlanineConfiguration(TMatConfiguration):
    def __init__(self):
        super().__init__()
        self.simulator = AlanineSimulator()


