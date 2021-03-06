"""Deal with parameters."""

import logging
import pickle
import numpy as np

log = logging.getLogger(__name__)


class AdaptiveParams:
    def __init__(self, spt, tpr, adapt_lt=0, build_lt=0, post_converge=0,
                 run_id=0, step_res=-1):
        self.tpr = tpr
        self.spt = spt
        self.run_id = run_id

        self.adapt_lt = adapt_lt
        self.build_lt = build_lt
        self.post_converge = post_converge

        self.step_res = step_res

    def save(self, fn):
        fn = "{}.pickl".format(fn)
        with open(fn, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, fn):
        with open(fn, 'rb') as f:
            return pickle.load(f)

    @property
    def dirname(self):
        return "blt-{build_lt}_alt-{adapt_lt}_spt-{spt}_tpr-{tpr}".format(
            build_lt=self.build_lt, adapt_lt=self.adapt_lt, spt=self.spt,
            tpr=self.tpr)

    @property
    def pretty_desc(self):
        """A description e.g. for plot titles."""
        return self.dirname

    @property
    def paramfn(self):
        return "params"

    @property
    def subbuild_uptos(self):
        if self.step_res <= 0:
            return np.asarray([self.spt])
        n_builds = self.spt // self.step_res
        assert self.spt % self.step_res == 0, 'Sub-build must divide'
        return (np.arange(n_builds) + 1) * self.step_res

