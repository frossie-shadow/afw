from __future__ import absolute_import, division, print_function

__all__ = []

from lsst.utils import continueClass
from ._fits import Fits


@continueClass
class Fits:
    def __enter__(self):
        return self

    def __exit__(self, cls, exc, traceback):
        self.closeFile()
