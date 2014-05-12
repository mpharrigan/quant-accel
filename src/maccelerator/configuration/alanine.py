

from maccelerator.configuration import base_classes as bc
import maccelerator as maccel


class AlanineSimulator(maccel.TMatSimulator):

    def __init__(self):
        # TODO: Pass in arguments
        super(maccel.TMatSimulator, self).__init__()


class AlanineConfiguration(bc.TMatConfiguration):

    def __init__(self):
        super(bc.TMatConfiguration, self).__init__()
        self.simulator = AlanineSimulator()


