"""A configuration that doesn't actually do anything,
but can be used for testing.
"""

import numpy as np

from ..simulate import Simulator
from ..model import Modeller, Model
from ..adapt import Adapter, SStates
from ..convergence.base import ConvergenceChecker, Convergence
from .base import Configuration
from ..param import AdaptiveParams


class SimpleSimulator(Simulator):
    """A simple simulation for testing."""

    def __init__(self):
        super().__init__()

    def simulate(self, sstate, params, traj_out_fn):
        """Output a list of numbers.

        :param sstate: Starting state (int)
        :param params: Contains number of steps to take
        :param traj_out_fn: Where to save the result

        """
        traj = np.arange(sstate, sstate + params.spt)

        if traj_out_fn is not None:
            np.save(traj_out_fn, traj)

    @property
    def trajfn(self):
        """Format string for where to save the trajectory."""
        return "traj-{traj_i}.npy"


class SimpleModel(Model):
    def __init__(self, trajs):
        super().__init__(msm=None)
        self.trajs = trajs


class SimpleModeller(Modeller):
    """A simple modeller for testing."""

    def __init__(self):
        super().__init__()
        self.trajs = []

    def model(self, traj_fns, params):
        """Just load up all of the 'trajectories'.

        :param params: Not used
        :param traj_fns: Files to load
        """
        trajs = [np.load(tfn) for tfn in traj_fns]
        return SimpleModel(trajs)


class SimpleConvchecker(ConvergenceChecker):
    """A simple check for convergence."""

    def __init__(self):
        super().__init__(tolerance=40)

    def check_convergence(self, model, params):
        """If we have enough trajectories, we're converged.

        :param model: The model for which we check convergence
        :param params: Not used
        """
        self.errors_over_time += [len(model.trajs)]
        converged = len(model.trajs) >= self.tolerance
        return Convergence(converged, self.errors_over_time)

    @property
    def n_plots(self):
        """We don't use plots for this dummy configuration."""
        return 0


class SimpleAdapter(Adapter):
    """A simple adapter for testing."""

    def adapt(self, model, params):
        """Get biggest number in the last n trajectories.

        :param model: The model from which we adapt
        :param params: Contains number of new states to generate
        """
        n_tpr = params.tpr
        ret = np.array([np.max(t) for t in model.trajs[-n_tpr:]])
        assert len(ret) == n_tpr
        return SStates(ret)

    def seed_states(self, params):
        """Start from zero

        :param params: Contains number of seed states to generate.
        """
        return SStates([0] * params.tpr)


class SimpleParams(AdaptiveParams):
    """Parameters for testing."""

    def __init__(self, spt, tpr):
        super().__init__(spt, tpr, run_id=0)

    @property
    def post_converge(self):
        """Do 4 rounds after convergence."""
        return 4

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
        self.adapter = SimpleAdapter()
        self.convchecker = SimpleConvchecker()
        self.simulator = SimpleSimulator()

    def get_param_grid(self):
        """Do two parameter configurations."""
        return [SimpleParams(spt=10, tpr=10), SimpleParams(spt=20, tpr=10)]
