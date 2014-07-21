"""Deal with parameters."""

from os.path import join as pjoin
import os
import logging

log = logging.getLogger()


class AdaptiveParams(object):
    def __init__(self, spt, tpr, run_id=0):
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

    @property
    def dirname(self):
        return "blt-{build_lt}_alt-{adapt_lt}_spt-{spt}_tpr-{tpr}".format(
            build_lt=self.build_lt, adapt_lt=self.adapt_lt, spt=self.spt,
            tpr=self.tpr
        )

    @property
    def pretty_desc(self):
        """A description e.g. for plot titles."""
        return self.dirname

