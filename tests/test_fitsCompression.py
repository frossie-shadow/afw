#
# LSST Data Management System
# Copyright 2017 LSST Corporation.
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
from builtins import range

import unittest

import numpy as np

import lsst.utils
import lsst.daf.base
import lsst.afw.geom
import lsst.afw.image
import lsst.afw.fits
import lsst.utils.tests


class ImageScalingTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(123, 456), lsst.afw.geom.Extent2I(78, 90))

    def tearDown(self):
        del self.bbox

    def checkRange(self, ImageClass, bitpix):
        image = ImageClass(self.bbox)
        mask = lsst.afw.image.Mask(self.bbox)
        bad = mask.getPlaneBitMask("BAD")
        base = 123
        swing = 100
        image.set(base)
        image.set(1, 1, base + swing)
        image.set(3, 3, 12345)
        mask.set(3, 3, bad)
        scaling = lsst.afw.fits.ImageScalingOptions(lsst.afw.fits.ImageScalingOptions.RANGE,
                                                    bitpix, [u"BAD"])
        # scale = scaling.determine(image, mask)
        # self.assertEqual(scale.bitpix, bitpix)
        # self.assertEqual(scale.bscale, swing/2**bitpix)

        # suffix = {16: "S", 32: "I", 64: "L"}[bitpix]
        # print(suffix)
        # disk = getattr(scale, "toDisk" + suffix)(image.getArray(), False)
        # expectType = {16: np.int16, 32: np.int32, 64: np.int64}[bitpix]
        # self.assertEqual(disk.dtype, expectType)
        # self.assertEqual(disk.min(), -2**(bitpix - 1))
        # self.assertEqual(disk.max(), 2**(bitpix - 1) - 1)

        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            with lsst.afw.fits.Fits(tmpFile, "w") as fits:
                options = lsst.afw.fits.ImageWriteOptions(scaling)
                header = lsst.daf.base.PropertyList()
                image.writeFits(fits, options, header)
            unpersisted = ImageClass(tmpFile)
            self.assertImagesAlmostEqual(image, unpersisted, rtol=2.0/2**bitpix)

    def checkStdev(self):
        pass


    def testRange(self):
        for cls in (lsst.afw.image.ImageF, lsst.afw.image.ImageD):
            for bitpix in (16, 32, 64):
                self.checkRange(cls, bitpix)


                    # lsst.afw.image.ImageU,
                    # lsst.afw.image.ImageI,
                    # lsst.afw.image.ImageL,
                    # lsst.afw.image.ImageF,
                    # lsst.afw.image.ImageD,

class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
