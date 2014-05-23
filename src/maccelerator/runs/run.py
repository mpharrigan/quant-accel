import os
import logging as log
from os.path import join as pjoin
from collections import defaultdict

from IPython.parallel import Client


TRAJ = 'trajs'
SSTATE = 'sstates'
MSMS = 'msms'
TRAJROUND = 'round-{round_i}'


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
        self.rundir = rundir
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

        if self.lbv is None:
            return False

        # Set up directories
        rundir = self.rundir
        traj_dir = pjoin(rundir, TRAJ)
        sstate_dir = pjoin(rundir, SSTATE)
        msms_dir = pjoin(rundir, MSMS)
        try:
            os.mkdir(rundir)
            os.mkdir(traj_dir)
            os.mkdir(sstate_dir)
            os.mkdir(msms_dir)
        except OSError as e:
            log.warning('Skipping %s (%s)', rundir, e)
            return False

        # Initialize variables for the loop
        converged = False
        round_i = 0
        rounds_left = self.params.post_converge
        sstate = self.config.modeller.seed_state(self.params.tpr)

        while True:
            # Make directory for trajectories for this round
            trajround_dir = pjoin(traj_dir, TRAJROUND.format(round_i=round_i))
            trajround_dir = os.path.abspath(trajround_dir)
            os.mkdir(trajround_dir)

            # Do fancy IPython.parallel stuff to parallelize simulation
            traj_outs = [pjoin(trajround_dir, self.trajfn.format(traj_i=i)) for
                         i in range(self.params.tpr)]
            self.lbv.map(self.config.simulator.simulate, sstate,
                         [self.params.spt] * self.params.tpr, traj_outs)
            self.trajs[round_i] = traj_outs

            # Model with all trajctories from this round and previous
            trajs_till_now = le_than(self.trajs, round_i)
            self.config.modeller.model(trajs_till_now)

            # Check convergence
            if not converged:
                converged = self.config.modeller.check_convergence()

            # Keep track of progress
            if converged:
                rounds_left -= 1

            # See if we're done
            if rounds_left <= 0:
                break

            # Move on
            sstate = self.config.adapter.adapt(self.params.tpr)
            round_i += 1


def le_than(traj_dict, round_i):
    """Pick trajectories whose round is less than or equal to round_i."""

    all_trajs = []
    for ri in traj_dict.keys():
        if ri <= round_i:
            all_trajs += traj_dict[ri]
    return all_trajs

