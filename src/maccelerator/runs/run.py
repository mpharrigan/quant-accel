"""Code to handle the adaptive loop."""

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
        self.trajs = defaultdict(list)
        self.rundir = rundir

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
        traj_dir, exit_status = self.config.file.make_directories(self.rundir)
        if not exit_status:
            return False

        # Initialize variables for the loop
        round_i = 0
        rounds_left = self.params.post_converge
        sstate = self.config.modeller.seed_state(self.params,
                                                 self.get_sstate_fn(-1))

        while True:
            log.info("Doing round %d", round_i)

            # Make directory for trajectories for this round
            trajround_dir = self.config.file.make_trajround(round_i)

            # Do fancy IPython.parallel stuff to parallelize simulation
            traj_outs = self.get_traj_fns(trajround_dir)
            self.lbv.map(self.config.simulator.simulate, sstate,
                         [self.params] * self.params.tpr, traj_outs)
            self.trajs[round_i] = traj_outs

            # Model with all trajectories from this round and previous
            trajs_till_now = le_than(self.trajs, round_i)
            m_success = self.config.modeller.model(trajs_till_now, self.params)

            # Get a new starting state
            ssout = self.get_sstate_fn(round_i)
            sstate = self.config.adapter.adapt(self.params, ssout)

            # Check convergence
            plot_fn = self.get_plot_fn(round_i)
            if m_success:
                converged = self.config.convchecker.check_convergence(
                    self.params)
                self.config.convchecker.plot_and_save(self.params, sstate,
                                                      plot_fn)
            else:
                converged = False
                self.config.convchecker.fallback(self.params)
                self.config.convchecker.plot_and_save(self.params, sstate,
                                                      plot_fn,
                                                      fallback=True)

            # Keep track of progress
            # Note: if we dip in and out of convergence it doesn't decrement
            # this, but it doesn't get reset.
            if converged:
                rounds_left -= 1

            # See if we're done
            if rounds_left < 0:
                break

            # Move on
            round_i += 1

    def get_traj_fns(self, trajround_dir):
        """Return a list of trajectory filenames for simulation."""
        return [
            pjoin(trajround_dir, self.config.simulator.trajfn.format(traj_i=i))
            for i in range(self.params.tpr)
        ]

    def get_sstate_fn(self, round_i):
        """Return a starting-state filename."""
        return pjoin(self.config.file.sstate_dir,
                     self.config.adapter.sstatefn.format(round_i=round_i + 1))

    def get_plot_fn(self, round_i):
        return pjoin(self.config.file.figs_dir,
                     self.config.convchecker.plotfn.format(round_i=round_i))


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
