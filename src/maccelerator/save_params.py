__author__ = 'harrigan'

import argparse
import json
import logging
import os

log = logging.getLogger()
log.setLevel(logging.INFO)


class ParamParser(object):

    def __init__(self, base_dir):
        self.base_dir = os.path.abspath(base_dir)
        self.param_locs = []

    def write_locations(self):
        """Save a list of all the locations of the param files.

        If this information isn't already stored in the object, generate it.
        """

        if len(self.param_locs) == 0:
            self.find_all()

        #TODO: write out a file

    def find_all(self):
        """Search for param files."""
        #TODO: Implement with subclasses probably

    def get_params_from_fs(self):
        """From the filesystem, get params and save them as json."""
        raise NotImplementedError("Subclass should implement this function")



class MullerParamParser(ParamParser):

    def get_params_from_fs(self):
        for dirpath, dirnames, filenames in os.walk(self.base_dir):
            log.debug("Looking in %s", dirpath)
            if ['msms', 'sstates', 'trajs'] in dirnames:
                log.info("Found a run directory %s", dirpath)

                cwd = os.path.abspath(dirpath)
                dnames = cwd.split('/')
                dname = dnames[-1]
                param_strs = dname.split('_')
                params = dict()
                for param_str in param_strs:
                    splits = param_str.split('-')
                    params[splits[0]] = splits[1]

                log.info("Params: %s", str(params))


def muller4_f(args):
    pp = MullerParamParser(args.base_dir)
    pp.get_params_from_fs()

def parse():
    parser = argparse.ArgumentParser(description="Save params as json",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--base_dir',
                        help='Directory to start looking from',
                        default='.')

    sp = parser.add_subparsers()
    muller4_p = sp.add_parser('muller4',
                              help="save parameters from file structure",
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    muller4_p.set_defaults(func=muller4_f)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    parse()
