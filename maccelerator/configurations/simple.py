"""A configuration that doesn't actually do anything,
but can be used for testing.
"""

import numpy as np

from ..simulate import Simulator
from ..model import Modeller, Model
from ..adapt import Adapter, SStates
from ..convergence.base import SubConvergenceChecker, SubConvergence, \
    SupConvergenceChecker
from .base import Configuration
from ..param import AdaptiveParams

from itertools import chain


class SimpleSimulator(Simulator):
    """A simple simulation for testing."""

    def __init__(self, config):
        super().__init__(config)

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
    """Just keep of the trajectories."""

    def __init__(self, params, trajs):
        super().__init__(params)
        self.trajs = trajs


class SimpleModeller(Modeller):
    """A simple modeller for testing."""

    def __init__(self, config):
        super().__init__(config)
        self.trajs = []

    def model(self, traj_fns, params, up_to=None):
        """Just load up all of the 'trajectories'.

        :param params: Not used
        :param traj_fns: Files to load
        :param up_to: Used for sub-modelling
        """
        trajs = [np.load(tfn) for tfn in chain(*traj_fns.values())]
        return SimpleModel(params, trajs)


class SimpleSupConvchecker(SupConvergenceChecker):
    """Wrapper for our one simple sub-convergence checker."""

    @classmethod
    def get_sub_checkers(cls, config):
        """Configure convergence criteria.
        :param config: Unused
        """
        return [SimpleSubConvchecker(tolerance=40)]


class SimpleSubConvchecker(SubConvergenceChecker):
    """A simple check for convergence."""

    def check_convergence(self, model, params):
        """If we have enough trajectories, we're converged.

        :param model: The model for which we check convergence
        :param params: Not used
        """
        self.errors_over_time += [len(model.trajs)]
        converged = len(model.trajs) >= self.tolerance
        return SubConvergence(converged, self.errors_over_time)

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
        super().__init__(spt, tpr, run_id=0, adapt_lt=1, build_lt=1,
                         post_converge=4)


class SimpleConfiguration(Configuration):
    """Encapsulate everything into a configuration object."""

    def defaults(config):
        """Set defaults."""

        # Simulator
        config.simulator_class = SimpleSimulator

        # Modeller
        config.modeller_class = SimpleModeller

        # Convergence checker
        config.convchecker_class = SimpleSupConvchecker

        # Adapter
        config.adapter_class = SimpleAdapter

        # Define a function to yield combinations of parameters
        def get_param_grid(run_id=0):
            """Return two sets of parameters."""
            return [SimpleParams(spt=10, tpr=10), SimpleParams(spt=20, tpr=10)]

        config.get_param_grid = get_param_grid

