"""Code to handle the adaptive loop."""

import logging
import pickle

from IPython.parallel import Client


log = logging.getLogger(__name__)


class MAccelRun(object):
    """Object for performing an accelerated run

    :param configuration: All the system-specific configuration
    :param params: The parameters (which are being varied across runs)
        for this specific run.
    :param rundir: Directory for this specific run.
    """

    def __init__(self, configuration, params, rundir, parallel=True):
        self.config = configuration
        self.params = params
        self.trajs = dict()
        self.rundir = rundir

        if parallel:
            try:
                c = Client()
                self.lbv = c.load_balanced_view()
                self.lbv.block = True
            except FileNotFoundError as e:
                log.error("Could not connect to parallel engine: %s", e)
                self.lbv = None
        else:
            self.lbv = NoParallelView()

    def run(self):
        """Execute adaptive loop until convergence."""

        # Aliases
        params = self.params
        config = self.config
        file = self.config.file
        config.convchecker.reset()

        # Check IPython.parallel
        if self.lbv is None:
            return False

        # Make directories
        exit_status = file.make_directories(self.rundir)
        if not exit_status:
            return False

        # Save parameters
        params.save(file.param_fn(params))

        # Initialize variables for the loop
        round_i = 0
        rounds_left = params.post_converge
        sstate = config.seed_state(params)
        sstate.save(file.sstate_fn(-1))

        while True:

            # Do fancy IPython.parallel stuff to parallelize simulation
            traj_is = range(params.tpr)
            traj_outs = file.make_traj_fns(round_i, traj_is)
            self.lbv.map(
                config.simulate,
                sstate.items(),
                [params] * params.tpr,
                traj_outs
            )
            self.trajs[round_i] = traj_outs

            # Make lots of models and check lots of convergence
            models = []
            converges = []
            for sub_i, up_to in enumerate(params.subbuild_uptos):
                model = config.model(self.trajs, params, up_to=up_to)
                converge = config.check_convergence(model, params)

                models.append(model)
                converges.append(converge)

            model = models[-1]
            converge = converges[-1]

            model_fn = "{}.pickl".format(file.model_fn(round_i))
            with open(model_fn, 'wb') as model_f:
                pickle.dump(models, model_f)

            conv_fn = "{}.pickl".format(file.conv_fn(round_i))
            with open(conv_fn, 'wb') as conv_f:
                pickle.dump(converges, conv_f)

            # Get a new starting state
            sstate = config.adapt(model, params)
            sstate.save(file.sstate_fn(round_i))

            log.info("Doing round %3d.\tConverged: %5s\tRounds left: %2d",
                     round_i, converge.converged, rounds_left)

            # Keep track of progress
            # Note: if we dip in and out of convergence it doesn't decrement
            # this, but it doesn't get reset.
            if converge.converged:
                rounds_left -= 1

            # See if we're done
            if rounds_left < 0:
                break

            # Move on
            round_i += 1

        # Save the run after done
        self.save(file.run_fn(self))


    def __getstate__(self):
        state = dict(self.__dict__)
        try:
            del state['lbv']
        except KeyError:
            pass
        return state

    @property
    def runfn(self):
        return 'run'

    def save(self, fn):
        fn = "{}.pickl".format(fn)
        with open(fn, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, fn):
        with open(fn, 'rb') as f:
            return pickle.load(f)

    @property
    def n_rounds(self):
        return len(self.trajs)


class NoParallelView:
    def __init__(self):
        pass

    def map(self, func, *args):
        for a1 in zip(*args):
            func(*a1)
