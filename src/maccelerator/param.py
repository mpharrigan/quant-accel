"""Deal with parameters."""


class AdaptiveParams(object):
    def __init__(self, spt, tpr, run_id):
        self.tpr = tpr
        self.spt = spt
        self.run_id = run_id

    @property
    def post_converge(self):
        raise NotImplementedError

    @property
    def adapt_lt(self):
        raise NotImplementedError

    @property
    def build_lt(self):
        raise NotImplementedError

