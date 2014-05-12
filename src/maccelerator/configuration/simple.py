"""A configuration that doesn't actually do anything,
but can be used for testing.
"""

from . import base_classes as bc
from ..simulate import Simulator
from ..model import Modeller


class SimpleSimulator(Simulator):
    def __init__(self):
        super().__init__()


class SimpleModeller(Modeller):
    def __init__(self):
        super().__init__()


class SimpleConfiguration(bc.Configuration):
    def __init__(self):
        super().__init__()
