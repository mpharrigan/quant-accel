__author__ = 'harrigan'

import argparse
import json
import logging
import os
import re

log = logging.getLogger()
log.addHandler(logging.StreamHandler())
log.setLevel(logging.INFO)


class ParamParser(object):
    def __init__(self, base_dir):
        self.base_dir = os.path.abspath(base_dir)
        self.param_locs = []

    def write_locations(self, loc_fn):
        """Save a list of all the locations of the param files.

        If this information isn't already stored in the object, generate it.
        """

        if len(self.param_locs) == 0:
            self.find_all()

        with open(loc_fn, 'w') as loc_f:
            loc_f.write('\n'.join(self.param_locs))


    def find_all(self):
        """Search for param files."""
        #TODO: Implement with subclasses probably
        raise NotImplementedError("Implement this later")

    def get_params_from_fs(self):
        """From the filesystem, get params and save them as json."""
        raise NotImplementedError("Subclass should implement this function")


class MullerParamParser(ParamParser):
    def get_params_from_fs(self, json_outfn):
        for dirpath, dirnames, filenames in os.walk(self.base_dir):
            log.debug("Looking in %s", dirpath)
            if 'msms' in dirnames and {'msms', 'sstates', 'trajs'} <= set(
                    dirnames):
                log.info("Found a run directory %s", dirpath)

                # Get params
                dnames = dirpath.split('/')
                dname = dnames[-1]
                param_strs = dname.split('_')
                params = dict()
                for param_str in param_strs:
                    splits = param_str.split('-')
                    params[splits[0]] = splits[1]

                match = re.search(r"runcopy-([0-9]+)", dirpath)

                # Try to figure out the runcopy number
                if match is not None:
                    params['runcopy'] = int(match.group(1))

                # Make everything into ints
                for par in params:
                    try:
                        params[par] = int(params[par])
                    except ValueError:
                        pass

                log.info("Params: %s", str(params))

                out_fn = os.path.join(dirpath, json_outfn)
                with open(out_fn, 'w') as json_f:
                    json.dump(params, json_f)

                self.param_locs += [out_fn]


def muller4_f(args):
    pp = MullerParamParser(args.base_dir)
    pp.get_params_from_fs('params.json')
    pp.write_locations('param_locs.dat')


def parse():
    parser = argparse.ArgumentParser(description="Save params as json",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('base_dir',
                        help='Directory to start looking from',
                        default='.')
    parser.add_argument('--debug', action='store_true')
    parser.set_defaults(debug=False)

    sp = parser.add_subparsers()
    muller4_p = sp.add_parser('muller4',
                              help="save parameters from file structure",
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    muller4_p.set_defaults(func=muller4_f)

    args = parser.parse_args()

    if args.debug:
        log.setLevel(logging.DEBUG)

    log.info("Info")
    log.debug("Debug")

    args.func(args)


if __name__ == "__main__":
    parse()
