'''
Created on Nov 11, 2013

@author: harrigan
'''
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

lagtimes = np.array(range(size))

my_results = list()

for i in range(2):
    my_results.append(("Hi %d" % rank, np.ones(4)*i + rank))


tots = None

comm.gather(my_results, root=0)

if rank == 0:
    print my_results
