__author__ = 'harrigan'

import maccelerator as maccel
from maccelerator.testing.utils import get_fn, get_folder
import os
import logging

logging.basicConfig(level=logging.DEBUG)

# Change directory
os.chdir('/home/harrigan/implement/quant-accel/july-tests/')

# Do it
configuration = maccel.AlanineConfiguration(get_fn('ala.msm.pickl'),
                                            get_fn('ala.centers.h5'))
param = maccel.AlanineParams(spt=200, tpr=1)
run = maccel.MAccelRun(configuration, param, get_folder('ala'))
run.run()


