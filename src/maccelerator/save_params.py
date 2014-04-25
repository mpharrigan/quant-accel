"""Previously, parameters (params) of simulations were strewn about, implicit
in directory names or otherwise. This file is to assist to moving to a scheme
where parameters are saved in a json file per run so there is no ambiguity."""

__author__ = 'harrigan'

import argparse
import json
import logging
import os
import re
import glob

log = logging.getLogger()
log.addHandler(logging.StreamHandler())
log.setLevel(logging.INFO)


class ParamParser(object):
    """Represents a parser that can take a directory structure and ouput
    param files (on a per-run basis) as well as list all the param files
    (and therefore all runs.)"""

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
    """Parse based on the directory structure used for the muller runs,
    specifically muller4."""


    # Deprecated
    def get_params_from_list(self, filelist_fn, json_outfn):
        """Run something like `find` to get a list of run dirs."""

        # The following is a valid command
        # However, globbing seems to be sufficient
        command = r"find -maxdetpth 3 -type d -name msms -exec dirname {} \; > rundirs.dat"

        with open(filelist_fn) as filelist_f:
            for line in filelist_f:
                line = line.strip()
                log.debug("Parsing directory %s", line)
                params = self.get_params_from_folder(line)
                out_fn = os.path.abspath(os.path.join(line, json_outfn))
                with open(out_fn, 'w') as json_f:
                    json.dump(params, json_f)
                self.param_locs += [out_fn]

    def get_params_from_fs(self, json_outfn):
        """From the filesystem, get params and save them as json.

        We expect a path {basedir}/runcopy-00/lt-20_spt-000_tpr-000/
        containing folder or file named msm

        This should be able to parse additional params if they are included
        in the folder name in a key-value_key2-value2 scheme.
        """

        dirglob = os.path.join(self.base_dir, "runcopy-*/*/msms")
        for dirpath in glob.iglob(dirglob):
            # Get rid of the trailing "msm"
            dirpath = os.path.dirname(dirpath)

            log.info("Found a run directory %s", dirpath)

            params = self.get_params_from_folder(dirpath)
            out_fn = os.path.abspath(os.path.join(dirpath, json_outfn))
            with open(out_fn, 'w') as json_f:
                json.dump(params, json_f)
            self.param_locs += [out_fn]


    def get_params_from_folder(self, dirpath):
        """Get params from a particular folder.

        We assume something else is doing a walk of the directory structure.
        """

        dirpath = os.path.abspath(dirpath)

        # Get params from folder name
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

        # Make everything into integers if possible
        for par in params:
            try:
                params[par] = int(params[par])
            except ValueError:
                pass

        log.info("Params: %s", str(params))
        return params


def muller4_f(args):
    """Entry point for argparse for doing muller work."""
    pp = MullerParamParser(args.base_dir)
    pp.get_params_from_fs('params.json')
    pp.write_locations('param_locs.dat')


def parse():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Save params as json",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('base_dir',
                        help='Directory to start looking from')
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
