"""Code to handle the adaptive loop."""

import os
import logging as log
from os.path import join as pjoin
from collections import defaultdict

from IPython.parallel import Client

log.basicConfig(level=log.DEBUG)


class MAccelRun(object):
    """Object for performing an accelerated run

    :param configuration: All the system-specific configuration
    :param params: The parameters (which are being varied across runs)
        for this specific run.
    :param rundir: Directory for this specific run.
    """

    def __init__(self, configuration, params, rundir):
        self.config = configuration
        self.params = params
        self.params.rundir = rundir
        self.trajs = defaultdict(list)

        # Template for trajectory names
        self.trajfn = configuration.simulator.trajfn

        try:
            c = Client()
            self.lbv = c.load_balanced_view()
            self.lbv.block = True
        except FileNotFoundError as e:
            log.error("Could not connect to parallel engine: %s", e)
            self.lbv = None

    def run(self):
        """Execute adaptive loop until convergence."""

        # Check IPython.parallel
        if self.lbv is None:
            return False

        # Make directories
        traj_dir, exit_status = self.params.make_directories()
        if not exit_status:
            return False

        # Initialize variables for the loop
        converged = False
        round_i = 0
        rounds_left = self.params.post_converge
        sstate = self.config.modeller.seed_state(self.params)

        while True:
            log.info("Doing round %d", round_i)

            # Make directory for trajectories for this round
            trajround_dir = self.params.make_trajround(round_i)

            # Do fancy IPython.parallel stuff to parallelize simulation
            traj_outs = [pjoin(trajround_dir, self.trajfn.format(traj_i=i)) for
                         i in range(self.params.tpr)]
            self.lbv.map(self.config.simulator.simulate, sstate,
                         [self.params] * self.params.tpr, traj_outs)
            self.trajs[round_i] = traj_outs

            # Model with all trajectories from this round and previous
            trajs_till_now = le_than(self.trajs, round_i)
            self.config.modeller.model(trajs_till_now, self.params)

            # Get a new starting state
            sstate = self.config.adapter.adapt(self.params)

            # Check convergence
            converged = self.config.convchecker.check_convergence(
                self.params, sstate)

            # Keep track of progress
            # Note: if we dip in and out of convergence it doesn't decrement
            # this, but it doesn't get reset.
            if converged:
                rounds_left -= 1

            # See if we're done
            if rounds_left <= 0:
                break

            # Move on
            round_i += 1


def le_than(traj_dict, round_i):
    """Pick trajectories whose round is less than or equal to round_i.

    :param traj_dict: key: round index (int); value: list of file names
    :param round_i: Index of round to get less than or equal to
    """

    all_trajs = []
    for ri in traj_dict.keys():
        if ri <= round_i:
            all_trajs += traj_dict[ri]
    return all_trajs
