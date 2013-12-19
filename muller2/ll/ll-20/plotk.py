import numpy as np
from matplotlib import pyplot as pp

ks = np.loadtxt('numberofclusters.dat')
pp.plot(ks[:,0], ks[:,1])

