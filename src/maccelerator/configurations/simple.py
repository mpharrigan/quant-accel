"""A configuration that doesn't actually do anything,
but can be used for testing.
"""

import numpy as np

from ..simulate import Simulator
from ..model import Modeller, Adapter
from ..configuration import Configuration
from ..param import AdaptiveParams


class SimpleSimulator(Simulator):
    """A simple simulation for testing."""

    def __init__(self):
        super().__init__()

    def simulate(self, sstate, n_steps, traj_out_fn):
        """Output a list of numbers.

        :param sstate: Starting state (int)
        :param n_steps: Number of steps to simulate
        :param traj_out_fn: Where to save the result

        :returns: A list or something
        """
        traj = range(sstate, n_steps)

        if traj_out_fn is not None:
            np.save(traj_out_fn, traj)
        return traj

    @property
    def trajfn(self):
        """Format string for where to save the trajectory."""
        return "traj-{traj_i}.npy"


class SimpleModeller(Modeller):
    """A simple modeller for testing."""

    def __init__(self):
        super().__init__()
        self.trajs = []

    def model(self, traj_fns):
        """Just load up all of the 'trajectories'.

        :param traj_fns: Files to load
        """
        self.trajs = [np.load(tfn) for tfn in traj_fns]

    def check_convergence(self):
        """If we have enough trajectories, we're converged."""
        return len(self.trajs) > 40

    def seed_state(self, tpr):
        """Start from 0 in each trajectory.

        :param tpr: Number of seed states to generate
        """
        return [0] * tpr


class SimpleAdapter(Adapter):
    """A simple adapter for testing."""

    def __init__(self, modeller):
        super().__init__()
        self.modeller = modeller

    def adapt(self, n_tpr):
        """Get biggest number in the last n trajectories.

        :param n_tpr: Number of new states to generate
        """
        ret = np.array([np.max(t) for t in self.modeller.trajs[-n_tpr:]])
        assert len(ret) == n_tpr
        return ret


class SimpleParams(AdaptiveParams):
    """Parameters for testing."""

    def __init__(self, spt, tpr):
        super().__init__(spt, tpr, run_id=0)

    @property
    def post_converge(self):
        """Do 5 rounds after convergence."""
        return 5

    @property
    def adapt_lt(self):
        """Dummy value."""
        return 1

    @property
    def build_lt(self):
        """Dummy value."""
        return 1


class SimpleConfiguration(Configuration):
    """Encapsulate all of the simple components for testing."""

    def __init__(self):
        super().__init__()

        self.modeller = SimpleModeller()
        self.adapter = SimpleAdapter(self.modeller)
        self.simulator = SimpleSimulator()

    def get_param_grid(self):
        """Do two parameter configurations."""
        return [
            SimpleParams(spt=10, tpr=10),
            SimpleParams(spt=20, tpr=10)
        ]
