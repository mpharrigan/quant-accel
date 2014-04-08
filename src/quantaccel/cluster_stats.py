'''
Created on Mar 20, 2014

@author: harrigan
'''

import argparse
from maccelerator import maccel
import numpy as np
import os
import mdtraj as md


def stats_for_round():
    _, trajs = maccel.load_trajectories(round_i)
    centroids = np.loadtxt(os.path.join('msms', 'centroids-%s-mk6-%d.npy' % (how, round_i)))
    
    

def parse():
    parser = argparse.ArgumentParser(help='Get stats for clustering')
    

if __name__ == "__main__":
    parse()
