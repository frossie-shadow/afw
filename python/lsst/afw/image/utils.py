#
# LSST Data Management System
# Copyright 2008-2017 LSST/AURA.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
from __future__ import absolute_import, division, print_function

__all__ = ["clipImage", "resetFilters", "defineFilter",
           "defineFiltersFromPolicy", "CalibNoThrow"]

from past.builtins import basestring
from builtins import object

import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDetect
from .maskedImage import MaskedImage, makeMaskedImage
from .mask import Mask
from .tanWcs import TanWcs
from .distortedTanWcs import DistortedTanWcs
from .filter import Filter, FilterProperty
from .calib import Calib


def clipImage(im, minClip, maxClip):
    """Clip an image to lie between minClip and maxclip (None to ignore)"""
    if isinstance(im, MaskedImage):
        mi = im
    else:
        mi = makeMaskedImage(im, Mask(im.getDimensions()))

    if minClip is not None:
        ds = afwDetect.FootprintSet(
            mi, afwDetect.Threshold(-minClip, afwDetect.Threshold.VALUE, False))
        afwDetect.setImageFromFootprintList(
            mi.getImage(), ds.getFootprints(), minClip)

    if maxClip is not None:
        ds = afwDetect.FootprintSet(mi, afwDetect.Threshold(maxClip))
        afwDetect.setImageFromFootprintList(
            mi.getImage(), ds.getFootprints(), maxClip)


def resetFilters():
    """Reset registry of filters and filter properties"""
    Filter.reset()
    FilterProperty.reset()


def defineFilter(name, lambdaEff, alias=[], force=False):
    """Define a filter and its properties in the filter registry"""
    prop = FilterProperty(name, lambdaEff, force)
    Filter.define(prop)
    if isinstance(alias, basestring):
        Filter.defineAlias(name, alias)
    else:
        for a in alias:
            Filter.defineAlias(name, a)


def defineFiltersFromPolicy(filterPolicy, reset=False):
    """Process a Policy and define the filters"""

    if reset:
        Filter.reset()
        FilterProperty.reset()
    #
    # Process the Policy and define the filters
    #
    policyFile = pexPolicy.DefaultPolicyFile(
        "afw", "FilterDictionary.paf", "policy")
    defPolicy = pexPolicy.Policy.createPolicy(
        policyFile, policyFile.getRepositoryPath(), True)

    filterPolicy.mergeDefaults(defPolicy.getDictionary())

    for p in filterPolicy.getArray("Filter"):
        Filter.define(FilterProperty(p.get("name"), p))
        if p.exists("alias"):
            for a in p.getArray("alias"):
                Filter.defineAlias(p.get("name"), a)


class CalibNoThrow(object):
    """A class intended to be used with python's with statement, to return NaNs for negative fluxes
    instead of raising exceptions (exceptions may be raised for other purposes).

E.g.
     with CalibNoThrow():
         ax.plot([exposure.getCalib().getMagnitude(a) for a in candAmps], zGood[:,k], 'b+')
    """

    def __enter__(self):
        self._throwOnNegative = Calib.getThrowOnNegativeFlux()
        Calib.setThrowOnNegativeFlux(False)

    def __exit__(self, *args):
        Calib.setThrowOnNegativeFlux(self._throwOnNegative)
