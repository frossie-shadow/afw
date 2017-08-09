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

import os
import unittest
import itertools

import numpy as np

import lsst.utils
import lsst.daf.base
import lsst.afw.geom
import lsst.afw.image
import lsst.afw.fits
import lsst.utils.tests

from lsst.afw.fits import ImageScalingOptions, ImageCompressionOptions


class ImageScalingTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(123, 456), lsst.afw.geom.Extent2I(7, 8))
        self.base = 456
        self.highValue = 789
        self.lowValue = 123
        self.maskedValue = 123456
        self.highPixel = (1, 1)
        self.lowPixel = (2, 2)
        self.maskedPixel = (3, 3)
        self.badMask = "BAD"
        self.stdev = 5.0

    def tearDown(self):
        del self.bbox

    def makeImage(self, ImageClass, scaling, addNoise=True):
        image = ImageClass(self.bbox)
        mask = lsst.afw.image.Mask(self.bbox)
        bad = mask.getPlaneBitMask(self.badMask)
        image.set(self.base)
        image.set(self.highPixel[0], self.highPixel[1], self.highValue)
        image.set(self.lowPixel[0], self.lowPixel[1], self.lowValue)
        image.set(self.maskedPixel[0], self.maskedPixel[1], self.maskedValue)
        mask.set(self.maskedPixel[0], self.maskedPixel[1], bad)

        rng = np.random.RandomState(12345)
        dtype = image.getArray().dtype
        if addNoise:
            image.getArray()[:] += rng.normal(0.0, self.stdev, image.getArray().shape).astype(dtype)

        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            with lsst.afw.fits.Fits(filename, "w") as fits:
                options = lsst.afw.fits.ImageWriteOptions(scaling)
                header = lsst.daf.base.PropertyList()
                image.writeFits(fits, options, header, mask)
            unpersisted = ImageClass(filename)
            self.assertEqual(image.getBBox(), unpersisted.getBBox())

            header = lsst.afw.image.readMetadata(filename)
            bscale = header.get("BSCALE")
            bzero = header.get("BZERO")

            self.assertEqual(header.get("BITPIX"), scaling.bitpix)

            if scaling.bitpix == 8:  # unsigned, says FITS
                maxValue = bscale*(2**scaling.bitpix) - 1 + bzero
                minValue = bzero
            else:
                maxValue = bscale*(2**(scaling.bitpix - 1) - 1) + bzero
                minValue = -bscale*2**(scaling.bitpix - 1) + bzero

        return image, unpersisted, bscale, bzero, minValue, maxValue

    def checkPixel(self, unpersisted, original, xy, expected, rtol=None, atol=None):
        self.assertFloatsAlmostEqual(unpersisted.get(*xy), expected, rtol=rtol, atol=atol)
        unpersisted.set(xy[0], xy[1], original.get(*xy))  # for ease of comparison of the whole image

    def checkSpecialPixels(self, original, unpersisted, maxValue, minValue, rtol=None, atol=None):
        highValue = original.get(*self.highPixel)
        lowValue = original.get(*self.lowPixel)
        maskedValue = original.get(*self.maskedPixel)

        expectHigh = min(highValue, maxValue)
        expectLow = max(lowValue, minValue)
        expectMasked = min(maskedValue, maxValue)

        self.checkPixel(unpersisted, original, self.highPixel, expectHigh, rtol=rtol, atol=atol)
        self.checkPixel(unpersisted, original, self.lowPixel, expectLow, rtol=rtol, atol=atol)
        self.checkPixel(unpersisted, original, self.maskedPixel, expectMasked, rtol=rtol, atol=atol)

    def checkRange(self, ImageClass, bitpix):
        print(ImageClass, bitpix)
        scaling = ImageScalingOptions(ImageScalingOptions.RANGE, bitpix, [u"BAD"], fuzz=False)
        original, unpersisted, bscale, bzero, minValue, maxValue = self.makeImage(ImageClass, scaling, False)

        bscaleExpect = (self.highValue - self.lowValue)/2**bitpix
        self.assertFloatsAlmostEqual(bscale, bscaleExpect, rtol=1.0e-6)

        rtol = 1.0/2**(bitpix - 1)
        self.checkSpecialPixels(original, unpersisted, maxValue, minValue, rtol=rtol)
        self.assertImagesAlmostEqual(original, unpersisted, rtol=rtol)

    def checkStdev(self, ImageClass, bitpix, scheme, addNoise=True):
        print(ImageClass, bitpix, scheme)

        quantizeLevel = 100.0
        quantizePad = 5.0
        scaling = lsst.afw.fits.ImageScalingOptions(scheme, bitpix, [u"BAD"], fuzz=False,
                                                    quantizeLevel=quantizeLevel, quantizePad=quantizePad)

        makeImageResults = self.makeImage(ImageClass, scaling)
        original, unpersisted, bscale, bzero, minValue, maxValue = makeImageResults

        self.assertFloatsAlmostEqual(bscale, self.stdev/quantizeLevel, rtol=3.0/quantizeLevel)

        # Testing bzero gives us the desired distribution
        if scheme == ImageScalingOptions.STDEV_POSITIVE:
            self.assertGreater(maxValue - self.base, self.base - minValue)
        elif scheme == ImageScalingOptions.STDEV_NEGATIVE:
            self.assertLess(maxValue - self.base, self.base - minValue)
        elif scheme == ImageScalingOptions.STDEV_BOTH:
            self.assertFloatsAlmostEqual(maxValue - self.base - 1, self.base - minValue, atol=1.0)
        else:
            raise RuntimeError("Unrecognised scheme: %d" % scheme)

        atol = 1.001*bscale
        self.checkSpecialPixels(original, unpersisted, maxValue, minValue, atol=atol)
        self.assertImagesAlmostEqual(original, unpersisted, atol=atol)

    def testRange(self):
        classList = (lsst.afw.image.ImageF, lsst.afw.image.ImageD)
        # NOT including BITPIX=64: provides larger dynamic range than 'double BSCALE' can handle
        bitpixList = (8, 16, 32)
        for cls, bitpix in itertools.product(classList, bitpixList):
            self.checkRange(cls, bitpix)

    def testStdev(self):
        classList = (lsst.afw.image.ImageF, lsst.afw.image.ImageD)
        # NOT including BITPIX=64: provides larger dynamic range than 'double BSCALE' can handle
        # NOT including BITPIX=8: tiny dynamic range means everything goes out of range easily
        bitpixList = (16, 32)
        schemeList = (ImageScalingOptions.STDEV_POSITIVE, ImageScalingOptions.STDEV_NEGATIVE,
                      ImageScalingOptions.STDEV_BOTH)
        for cls, bitpix, scheme in itertools.product(classList, bitpixList, schemeList):
            self.checkStdev(cls, bitpix, scheme)

    def testRangeFailures(self):
        classList = (lsst.afw.image.ImageU, lsst.afw.image.ImageI, lsst.afw.image.ImageL)
        bitpixList = (8, 16, 32)
        for cls, bitpix in itertools.product(classList, bitpixList):
            with self.assertRaises(lsst.pex.exceptions.InvalidParameterError):
                self.checkRange(cls, bitpix)

    def testStdevFailures(self):
        classList = (lsst.afw.image.ImageU, lsst.afw.image.ImageI, lsst.afw.image.ImageL)
        bitpixList = (16, 32)
        schemeList = (ImageScalingOptions.STDEV_POSITIVE, ImageScalingOptions.STDEV_NEGATIVE,
                      ImageScalingOptions.STDEV_BOTH)
        for cls, bitpix, scheme in itertools.product(classList, bitpixList, schemeList):
            with self.assertRaises(lsst.pex.exceptions.InvalidParameterError):
                # Adding noise to integer image results in:
                #     TypeError: Cannot cast ufunc add output from dtype('float64') to
                #     dtype('uint16') with casting rule 'same_kind'
                self.checkStdev(cls, bitpix, scheme, addNoise=False)


class ImageCompressionTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(123, 456), lsst.afw.geom.Extent2I(78, 90))
        self.background = 12345.6789
        self.noise = 67.89

    def makeImage(self, ImageClass):
        image = ImageClass(self.bbox)
        rng = np.random.RandomState(12345)
        dtype = image.getArray().dtype
        noise = rng.normal(0.0, self.noise, image.getArray().shape).astype(dtype)
        image.getArray()[:] = np.array(self.background, dtype=dtype) + noise
        return image

    def makeMask(self):
        mask = lsst.afw.image.Mask(self.bbox)
        rng = np.random.RandomState(12345)
        dtype = mask.getArray().dtype
        # Note: this is not the usual distribution of mask pixel values, so the compression factors
        # we get here will not be indicative of what we will get in practise.
        mask.getArray()[:] = rng.randint(0, 2**(dtype.itemsize*8 - 1), mask.getArray().shape, dtype=dtype)
        return mask

    def checkCompressedImage(self, ImageClass, image, compression, scaling=None, atol=0.0):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            with lsst.afw.fits.Fits(filename, "w") as fits:
                if scaling:
                    options = lsst.afw.fits.ImageWriteOptions(compression, scaling)
                else:
                    options = lsst.afw.fits.ImageWriteOptions(compression)
                header = lsst.daf.base.PropertyList()
                image.writeFits(fits, options, header)
                fileSize = os.stat(filename).st_size
                numBlocks = 1 + np.ceil(self.bbox.getArea()*image.getArray().dtype.itemsize/2880.0)
                uncompressedSize = 2880*numBlocks
                print(ImageClass, compression.scheme, fileSize, uncompressedSize, fileSize/uncompressedSize)

            unpersisted = ImageClass(filename)
            self.assertEqual(image.getBBox(), unpersisted.getBBox())
            self.assertImagesAlmostEqual(unpersisted, image, atol=atol)

    def testLosslessFloat(self):
        classList = (lsst.afw.image.ImageF, lsst.afw.image.ImageD)
        schemeList = ("GZIP", "GZIP_SHUFFLE")  # Lossless float compression requires GZIP
        for cls, scheme in itertools.product(classList, schemeList):
            image = self.makeImage(cls)
            compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme))
            self.checkCompressedImage(cls, image, compression, atol=0.0)

    def testLosslessInt(self):
        # Compression of ImageL is unsupported by cfitsio
        classList = (lsst.afw.image.ImageU, lsst.afw.image.ImageI)
        schemeList = ("GZIP", "GZIP_SHUFFLE", "RICE")
        for cls, scheme in itertools.product(classList, schemeList):
            if scheme == "RICE":
                # cfitsio 3.36 requres the quantizeLevel to be set, even for an integer image, to
                # avoid the following error:
                #
                # error compressing image (413)
                #   Lossless compression of floating point images must use GZIP (imcomp_init_table)
                compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme),
                                                      1.2345)
            else:
                compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme))

            image = self.makeImage(cls)
            self.checkCompressedImage(cls, image, compression, atol=0.0)

    def testMask(self):
        # Default mask type (32 bits) is out of range for PLIO (24 bits)
        for scheme in ("GZIP", "GZIP_SHUFFLE", "RICE"):
            if scheme == "RICE":
                # cfitsio 3.36 requres the quantizeLevel to be set, even for an integer image, to
                # avoid the following error:
                #
                # error compressing image (413)
                #   Lossless compression of floating point images must use GZIP (imcomp_init_table)
                compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme),
                                                      1.2345)
            else:
                compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme))

            self.checkCompressedImage(lsst.afw.image.Mask, self.makeMask(), compression, atol=0.0)

    def testLossyFloatCfitsio(self):
        # cfitsio lossy float compression is controlled through the quantize_level parameter
        classList = (lsst.afw.image.ImageF, lsst.afw.image.ImageD)
        schemeList = ("GZIP", "GZIP_SHUFFLE", "RICE")
        quantizeList = (4.0, 10.0)
        for cls, scheme, quantizeLevel in itertools.product(classList, schemeList, quantizeList):
            compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme),
                                                  quantizeLevel)
            image = self.makeImage(cls)
            self.checkCompressedImage(cls, image, compression, atol=1.25*self.noise/quantizeLevel)

    def testLossyFloatOurs(self):
        # Our lossy float compression is done by scaling first
        classList = (lsst.afw.image.ImageF, lsst.afw.image.ImageD)
        schemeList = ("GZIP", "GZIP_SHUFFLE", "RICE")
        bitpixList = (16, 32)
        quantizeList = (4.0, 10.0)
        for cls, scheme, bitpix, quantize, fuzz in itertools.product(classList, schemeList, bitpixList,
                                                                     quantizeList, [False, True]):
            cfitsioQuantize = 1.2345 if scheme == "RICE" else 0.0
            compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme),
                                                  cfitsioQuantize)
            scaling = ImageScalingOptions(ImageScalingOptions.STDEV_BOTH, bitpix, quantizeLevel=quantize,
                                          fuzz=fuzz)
            image = self.makeImage(cls)
            self.checkCompressedImage(cls, image, compression, scaling, atol=1.2*self.noise/quantize)

    def checkCompressedMaskedImage(self, image, imageOptions, maskOptions, varianceOptions, atol=0.0):
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            with lsst.afw.fits.Fits(filename, "w") as fits:
                image.writeFits(fits, imageOptions, maskOptions, varianceOptions)
            unpersisted = type(image)(filename)
            if hasattr(image, "getMaskedImage"):
                image = image.getMaskedImage()
                unpersisted = unpersisted.getMaskedImage()
            self.assertEqual(image.getBBox(), unpersisted.getBBox())
            self.assertImagesAlmostEqual(unpersisted.getImage(), image.getImage(), atol=atol)
            self.assertImagesAlmostEqual(unpersisted.getMask(), image.getMask(), atol=atol)
            self.assertImagesAlmostEqual(unpersisted.getVariance(), image.getVariance(), atol=atol)

    def checkMaskedImage(self, imageOptions, maskOptions, varianceOptions, atol=0.0):
        image = lsst.afw.image.makeMaskedImage(self.makeImage(lsst.afw.image.ImageF),
                                               self.makeMask(), self.makeImage(lsst.afw.image.ImageF))
        self.checkCompressedMaskedImage(image, imageOptions, maskOptions, varianceOptions, atol=atol)
        exp = lsst.afw.image.makeExposure(image)
        self.checkCompressedMaskedImage(exp, imageOptions, maskOptions, varianceOptions, atol=atol)

    def testMaskedImage(self):
        # Lossless
        lossless = lsst.afw.fits.ImageCompressionOptions(ImageCompressionOptions.GZIP_SHUFFLE)
        options = lsst.afw.fits.ImageWriteOptions(lossless)
        self.checkMaskedImage(options, options, options, atol=0.0)

        # Lossy cfitsio compression
        quantize = 4.0
        cfitsio = lsst.afw.fits.ImageCompressionOptions(ImageCompressionOptions.GZIP_SHUFFLE, quantize)
        imageOptions = lsst.afw.fits.ImageWriteOptions(cfitsio)
        maskOptions = lsst.afw.fits.ImageWriteOptions(lossless)
        self.checkMaskedImage(imageOptions, maskOptions, imageOptions, atol=1.25*self.noise/quantize)

        # Lossy our compression
        quantize = 4.0
        compression = lsst.afw.fits.ImageCompressionOptions(ImageCompressionOptions.GZIP_SHUFFLE, 0.0)
        scaling = lsst.afw.fits.ImageScalingOptions(ImageScalingOptions.STDEV_BOTH, 32,
                                                    quantizeLevel=quantize)
        imageOptions = lsst.afw.fits.ImageWriteOptions(compression, scaling)
        maskOptions = lsst.afw.fits.ImageWriteOptions(compression)
        self.checkMaskedImage(imageOptions, maskOptions, imageOptions, atol=1.2*self.noise/quantize)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    import sys
    setup_module(sys.modules[__name__])
    unittest.main()
