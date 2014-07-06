__author__ = 'harrigan'


import maccelerator as maccel
import os
import logging

logging.basicConfig(level=logging.DEBUG)

# Change directory
os.chdir('/home/harrigan/implement/quant-accel/july-tests/')
sys_xml_fn = '/home/harrigan/projects/quant-accel/reference/muller_sys.xml'
int_xml_fn = '/home/harrigan/projects/quant-accel/reference/muller_int.xml'

# Make a directory
i = 0
while True:
    try:
        folder_name = 'r{}'.format(i)
        os.mkdir(folder_name)
        break
    except:
        i+=1


configuration = maccel.MullerConfiguration(sys_xml_fn, int_xml_fn)
param = maccel.MullerParams(spt=800, tpr=2)
run = maccel.MAccelRun(configuration, param, folder_name)
run.run()

