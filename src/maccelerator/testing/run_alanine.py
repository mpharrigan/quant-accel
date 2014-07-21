__author__ = 'harrigan'

import maccelerator as maccel
import os
import logging

logging.basicConfig(level=logging.DEBUG)

# Change directory
os.chdir('/home/harrigan/implement/quant-accel/july-tests/')
alanine_refmsm_fn = '/home/harrigan/projects/quant-accel/reference/ala.msm.pickl'
alanine_centers_fn = '/home/harrigan/projects/quant-accel/reference/ala.centers.h5'

# Make a directory
i = 0
while True:
    try:
        folder_name = 'a{}'.format(i)
        os.mkdir(folder_name)
        break
    except:
        i += 1

configuration = maccel.AlanineConfiguration(alanine_refmsm_fn,
                                            alanine_centers_fn)
param = maccel.AlanineParams(spt=100, tpr=1)
run = maccel.MAccelRun(configuration, param, folder_name)
run.run()


