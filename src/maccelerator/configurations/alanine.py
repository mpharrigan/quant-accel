from ..simulate import TMatSimulator
from ..model import TMatModeller, SortCountsAdapter
from ..configuration import TMatConfiguration
from ..param import AdaptiveParams


class AlanineSimulator(TMatSimulator):
    pass


class AlanineModeller(TMatModeller):
    pass


class AlanineParams(AdaptiveParams):
    @property
    def build_lt(self):
        return 1

    @property
    def adapt_lt(self):
        return 1


class AlanineAdapter(SortCountsAdapter):
    pass


class AlanineConfiguration(TMatConfiguration):
    def __init__(self, tmat_fn):
        super().__init__(tmat_fn)
        self.simulator = AlanineSimulator(self.tmat)
        self.modeller = AlanineModeller()
        self.adapter = AlanineAdapter(self.modeller)


