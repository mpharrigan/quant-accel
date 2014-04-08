'''
Created on Mar 5, 2014

@author: harrigan
'''

import argparse
import os
from matplotlib import pyplot as pp
import mdtraj as md
from toy_accel import mullerforce as mf


def plot_trajs(traj_dir):
    trajs = [md.load(os.path.join(traj_dir, s))
             for s in os.listdir(traj_dir) if s.endswith('.h5')]
    mf.MullerForce.plot()
    for t in trajs:
        pp.plot(t.xyz[:, 0, 0], t.xyz[:, 0, 1], 'o-')
    pp.show()


def parse():
    parser = argparse.ArgumentParser(description='Perform accelerated sampling',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '--dir',
        '-d',
        help='directory to look for trajs',
        default='.')
    args = parser.parse_args()
    plot_trajs(args.dir)

if __name__ == "__main__":
    parse()
