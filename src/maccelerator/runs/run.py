

class MAccelRun(object):
    """Object for performing an accelerated run

    :param configuration: All the system-specific configuration
    :param params: The parameters (which are being varied across runs)
        for this specific run.
    """
    def __init__(self, configuration, params):
        self.config = configuration
        self.params = params

    def run(self):

        converged = False
        rounds_left = self.params.post_converge
        sstate = self.config.modeller.seed_state()

        while True:

            for i in range(self.params.tpr):
                self.config.simulator.simulate(sstate, self.params.spt)

            self.config.modeller.model()

            if not converged:
                converged = self.config.modeller.check_convergence()

            if converged:
                rounds_left -= 1

            if rounds_left <= 0:
                break

            sstate = self.config.modeller.adapt()

