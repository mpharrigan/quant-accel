from ..simulate import TMatSimulator
from ..configuration import TMatConfiguration


class AlanineSimulator(TMatSimulator):
    def __init__(self):
        # TODO: Pass in arguments
        super().__init__()


class AlanineConfiguration(TMatConfiguration):
    def __init__(self):
        super().__init__()
        self.simulator = AlanineSimulator()


