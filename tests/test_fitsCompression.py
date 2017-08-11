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
    """Tests of image scaling

    The pattern here is to create an image, write it out with a
    specific scaling scheme, read it back in and test that everything
    is as we expect. We do this for each scaling scheme in its own
    test, and within that test iterate over various parameters (input
    image type, BITPIX, etc.). The image we create has a few features
    (low, high and masked pixels) that we check.
    """
    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(123, 456), lsst.afw.geom.Extent2I(7, 8))
        self.base = 456  # Base value for pixels
        self.highValue = 789  # Value for high pixel
        self.lowValue = 123  # Value for low pixel
        self.maskedValue = 123456  # Value for masked pixel (to throw off statistics)
        self.highPixel = (1, 1)  # Location of high pixel
        self.lowPixel = (2, 2)  # Location of low pixel
        self.maskedPixel = (3, 3)  # Location of masked pixel
        self.badMask = "BAD"  # Mask plane to set for masked pixel
        self.stdev = 5.0  # Noise stdev to add to image

    def makeImage(self, ImageClass, scaling, addNoise=True):
        """Make an image for testing

        We create an image, persist and unpersist it, returning
        some data to the caller.

        Parameters
        ----------
        ImageClass : `type`, an `lsst.afw.image.Image` class
            Class of image to create.
        scaling : `lsst.afw.fits.ImageScalingOptions`
            Scaling to apply during persistence.
        addNoise : `bool`
            Add noise to image?

        Returns
        -------
        image : `lsst.afw.image.Image` (ImageClass)
            Created image.
        unpersisted : `lsst.afw.image.Image` (ImageClass)
            Unpersisted image.
        bscale, bzero : `float`
            FITS scale factor and zero used.
        minValue, maxValue : `float`
            Minimum and maximum value given the nominated scaling.
        """
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

            if scaling.scheme != ImageScalingOptions.NONE:
                self.assertEqual(header.get("BITPIX"), scaling.bitpix)

            if scaling.bitpix == 8:  # unsigned, says FITS
                maxValue = bscale*(2**scaling.bitpix) - 1 + bzero
                minValue = bzero
            else:
                maxValue = bscale*(2**(scaling.bitpix - 1) - 1) + bzero
                minValue = -bscale*2**(scaling.bitpix - 1) + bzero

        return image, unpersisted, bscale, bzero, minValue, maxValue

    def checkPixel(self, unpersisted, original, xy, expected, rtol=None, atol=None):
        """Check one of the special pixels

        After checking, we set this pixel to the original value so
        it's then easy to compare the entire image.

        Parameters
        ----------
        unpersisted : `lsst.afw.image.Image`
            Unpersisted image.
        original : `lsst.afw.image.Image`
            Original image.
        xy : `tuple` of two `int`s
            Position of pixel to check.
        expected : scalar
            Expected value of pixel.
        rtol, atol : `float` or `None`
            Relative/absolute tolerance for comparison.
        """
        self.assertFloatsAlmostEqual(unpersisted.get(*xy), expected, rtol=rtol, atol=atol)
        unpersisted.set(xy[0], xy[1], original.get(*xy))  # for ease of comparison of the whole image

    def checkSpecialPixels(self, original, unpersisted, maxValue, minValue, rtol=None, atol=None):
        """Check the special pixels

        Parameters
        ----------
        original : `lsst.afw.image.Image`
            Original image.
        unpersisted : `lsst.afw.image.Image`
            Unpersisted image.
        minValue, maxValue : `float`
            Minimum and maximum value given the nominated scaling.
        rtol, atol : `float` or `None`
            Relative/absolute tolerance for comparison.
        """
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
        """Check that the RANGE scaling works

        Parameters
        ----------
        ImageClass : `type`, an `lsst.afw.image.Image` class
            Class of image to create.
        bitpix : `int`
            Bits per pixel for FITS image.
        """
        scaling = ImageScalingOptions(ImageScalingOptions.RANGE, bitpix, [u"BAD"], fuzz=False)
        original, unpersisted, bscale, bzero, minValue, maxValue = self.makeImage(ImageClass, scaling, False)

        bscaleExpect = (self.highValue - self.lowValue)/2**bitpix
        self.assertFloatsAlmostEqual(bscale, bscaleExpect, rtol=1.0e-6)

        rtol = 1.0/2**(bitpix - 1)
        self.checkSpecialPixels(original, unpersisted, maxValue, minValue, rtol=rtol)
        self.assertImagesAlmostEqual(original, unpersisted, rtol=rtol)

    def checkStdev(self, ImageClass, bitpix, scheme):
        """Check that one of the STDEV scaling schemes work

        Parameters
        ----------
        ImageClass : `type`, an `lsst.afw.image.Image` class
            Class of image to create.
        bitpix : `int`
            Bits per pixel for FITS image.
        scheme : `lsst.afw.fits.ImageScalingOptions.ScalingScheme`
            Scaling scheme to apply (one of the STDEV_*).
        """
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
        """Test that the RANGE scaling works on floating-point inputs

        We deliberately don't include BITPIX=64 because int64 provides
        a larger dynamic range than 'double BSCALE' can handle.
        """
        classList = (lsst.afw.image.ImageF, lsst.afw.image.ImageD)
        bitpixList = (8, 16, 32)
        for cls, bitpix in itertools.product(classList, bitpixList):
            self.checkRange(cls, bitpix)

    def testStdev(self):
        """Test that the STDEV scalings work on floating-point inputs

        We deliberately don't include BITPIX=64 because int64 provides
        a larger dynamic range than 'double BSCALE' can handle.

        We deliberately don't include BITPIX=8 because that provides
        only a tiny dynamic range where everything goes out of range easily.
        """
        classList = (lsst.afw.image.ImageF, lsst.afw.image.ImageD)
        bitpixList = (16, 32)
        schemeList = (ImageScalingOptions.STDEV_POSITIVE, ImageScalingOptions.STDEV_NEGATIVE,
                      ImageScalingOptions.STDEV_BOTH)
        for cls, bitpix, scheme in itertools.product(classList, bitpixList, schemeList):
            self.checkStdev(cls, bitpix, scheme)

    def testRangeFailures(self):
        """Test that the RANGE scaling fails on integer inputs"""
        classList = (lsst.afw.image.ImageU, lsst.afw.image.ImageI, lsst.afw.image.ImageL)
        bitpixList = (8, 16, 32)
        for cls, bitpix in itertools.product(classList, bitpixList):
            with self.assertRaises(lsst.pex.exceptions.InvalidParameterError):
                self.checkRange(cls, bitpix)

    def testStdevFailures(self):
        """Test that the STDEV scalings fail on integer inputs"""
        classList = (lsst.afw.image.ImageU, lsst.afw.image.ImageI, lsst.afw.image.ImageL)
        bitpixList = (16, 32)
        schemeList = (ImageScalingOptions.STDEV_POSITIVE, ImageScalingOptions.STDEV_NEGATIVE,
                      ImageScalingOptions.STDEV_BOTH)
        for cls, bitpix, scheme in itertools.product(classList, bitpixList, schemeList):
            with self.assertRaises(lsst.pex.exceptions.InvalidParameterError):
                self.checkStdev(cls, bitpix, scheme)

    def checkNone(self, ImageClass, bitpix):
        """Check that the NONE scaling scheme works

        Parameters
        ----------
        ImageClass : `type`, an `lsst.afw.image.Image` class
            Class of image to create.
        bitpix : `int`
            Bits per pixel for FITS image.
        """
        scaling = ImageScalingOptions(ImageScalingOptions.NONE, bitpix, [u"BAD"], fuzz=False)
        original, unpersisted, bscale, bzero, minValue, maxValue = self.makeImage(ImageClass, scaling)
        self.assertFloatsAlmostEqual(bscale, 1.0, atol=0.0)
        self.assertFloatsAlmostEqual(bzero, 0.0, atol=0.0)
        self.assertImagesAlmostEqual(original, unpersisted, atol=0.0)

    def testNone(self):
        """Test that the NONE scaling works on floating-point inputs"""
        classList = (lsst.afw.image.ImageF, lsst.afw.image.ImageD)
        bitpixList = (8, 16, 32)
        for cls, bitpix in itertools.product(classList, bitpixList):
            self.checkNone(cls, bitpix)

    def checkManual(self, ImageClass, bitpix):
        """Check that the MANUAL scaling scheme works

        Parameters
        ----------
        ImageClass : `type`, an `lsst.afw.image.Image` class
            Class of image to create.
        bitpix : `int`
            Bits per pixel for FITS image.
        """
        bscaleSet = 1.2345
        bzeroSet = self.base
        scaling = ImageScalingOptions(ImageScalingOptions.MANUAL, bitpix, [u"BAD"], bscale=bscaleSet,
                                      bzero=bzeroSet, fuzz=False)
        original, unpersisted, bscale, bzero, minValue, maxValue = self.makeImage(ImageClass, scaling)
        self.assertFloatsAlmostEqual(bscale, bscaleSet, atol=0.0)
        self.assertFloatsAlmostEqual(bzero, bzeroSet, atol=0.0)
        self.assertImagesAlmostEqual(original, unpersisted, atol=bscale)

    def testManual(self):
        """Test that the MANUAL scaling works on floating-point inputs"""
        classList = (lsst.afw.image.ImageF, lsst.afw.image.ImageD)
        bitpixList = (16, 32)
        for cls, bitpix in itertools.product(classList, bitpixList):
            self.checkNone(cls, bitpix)


class ImageCompressionTestCase(lsst.utils.tests.TestCase):
    """Tests of image compression

    We test compression both with and without loss (quantisation/scaling).

    The pattern here is to create an image, write it out with a
    specific compression scheme, read it back in and test that everything
    is as we expect. We do this for each compression scheme in its own
    test, and within that test iterate over various parameters (input
    image type, BITPIX, etc.).

    We print the (inverse) compression ratio for interest. Note that
    these should not be considered to be representative of the
    compression that will be achieved on scientific data, since the
    images created here have different qualities than scientific data
    that will affect the compression ratio (e.g., size, noise properties).
    """
    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(123, 456), lsst.afw.geom.Extent2I(78, 90))
        self.background = 12345.6789  # Background value
        self.noise = 67.89  # Noise (stdev)

    def makeImage(self, ImageClass):
        """Create an image

        Parameters
        ----------
        ImageClass : `type`, an `lsst.afw.image.Image` class
            Class of image to create.

        Returns
        -------
        image : `ImageClass`
            The created image.
        """
        image = ImageClass(self.bbox)
        rng = np.random.RandomState(12345)
        dtype = image.getArray().dtype
        noise = rng.normal(0.0, self.noise, image.getArray().shape).astype(dtype)
        image.getArray()[:] = np.array(self.background, dtype=dtype) + noise
        return image

    def makeMask(self):
        """Create a mask

        Note that we generate a random distribution of mask pixel values,
        which is very different from the usual distribution in science images.

        Returns
        -------
        mask : `lsst.afw.image.Mask`
            The created mask.
        """
        mask = lsst.afw.image.Mask(self.bbox)
        rng = np.random.RandomState(12345)
        dtype = mask.getArray().dtype
        mask.getArray()[:] = rng.randint(0, 2**(dtype.itemsize*8 - 1), mask.getArray().shape, dtype=dtype)
        return mask

    def checkCompressedImage(self, ImageClass, image, compression, scaling=None, atol=0.0):
        """Check that compression works on an image

        Parameters
        ----------
        ImageClass : `type`, an `lsst.afw.image.Image` class
            Class of image.
        image : `lsst.afw.image.Image`
            Image to compress.
        compression : `lsst.afw.fits.ImageCompressionOptions`
            Compression parameters.
        scaling : `lsst.afw.fits.ImageScalingOptions` or `None`
            Scaling parameters for lossy compression (optional).
        atol : `float`
            Absolute tolerance for comparing unpersisted image.
        """
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
        """Test lossless compression of floating-point image"""
        classList = (lsst.afw.image.ImageF, lsst.afw.image.ImageD)
        schemeList = ("GZIP", "GZIP_SHUFFLE")  # Lossless float compression requires GZIP
        for cls, scheme in itertools.product(classList, schemeList):
            image = self.makeImage(cls)
            compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme))
            self.checkCompressedImage(cls, image, compression, atol=0.0)

    def testLosslessInt(self):
        """Test lossless compression of integer image

        We deliberately don't test `lsst.afw.image.ImageL` because
        compression of LONGLONG images is unsupported by cfitsio.
        """
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
                                                      quantizeLevel=1.2345)
            else:
                compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme))

            image = self.makeImage(cls)
            self.checkCompressedImage(cls, image, compression, atol=0.0)

    def testMask(self):
        """Test compression of mask

        We deliberately don't test PLIO compression (which is designed for
        masks) because our default mask type (32) has too much dynamic range
        for PLIO (limit of 24 bits).
        """
        for scheme in ("GZIP", "GZIP_SHUFFLE", "RICE"):
            if scheme == "RICE":
                # cfitsio 3.36 requres the quantizeLevel to be set, even for an integer image, to
                # avoid the following error:
                #
                # error compressing image (413)
                #   Lossless compression of floating point images must use GZIP (imcomp_init_table)
                compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme),
                                                      quantizeLevel=1.2345)
            else:
                compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme))

            self.checkCompressedImage(lsst.afw.image.Mask, self.makeMask(), compression, atol=0.0)

    def testLossyFloatCfitsio(self):
        """Test lossy compresion of floating-point images with cfitsio

        cfitsio does the compression, controlled through the 'quantizeLevel'
        parameter. Note that cfitsio doesn't have access to our masks when
        it does its statistics.
        """
        classList = (lsst.afw.image.ImageF, lsst.afw.image.ImageD)
        schemeList = ("GZIP", "GZIP_SHUFFLE", "RICE")
        quantizeList = (4.0, 10.0)
        for cls, scheme, quantizeLevel in itertools.product(classList, schemeList, quantizeList):
            compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme),
                                                  quantizeLevel=quantizeLevel)
            image = self.makeImage(cls)
            self.checkCompressedImage(cls, image, compression, atol=1.25*self.noise/quantizeLevel)

    def testLossyFloatOurs(self):
        """Test lossy compression of floating-point images ourselves

        We do lossy compression by scaling first. We have full control over
        the scaling (multiple scaling schemes), and we have access to our
        own masks when we do statistics.
        """
        classList = (lsst.afw.image.ImageF, lsst.afw.image.ImageD)
        schemeList = ("GZIP", "GZIP_SHUFFLE", "RICE")
        bitpixList = (16, 32)
        quantizeList = (4.0, 10.0)
        for cls, scheme, bitpix, quantize, fuzz in itertools.product(classList, schemeList, bitpixList,
                                                                     quantizeList, [False, True]):
            cfitsioQuantize = 1.2345 if scheme == "RICE" else 0.0
            compression = ImageCompressionOptions(lsst.afw.fits.compressionSchemeFromString(scheme),
                                                  quantizeLevel=cfitsioQuantize)
            scaling = ImageScalingOptions(ImageScalingOptions.STDEV_BOTH, bitpix, quantizeLevel=quantize,
                                          fuzz=fuzz)
            image = self.makeImage(cls)
            self.checkCompressedImage(cls, image, compression, scaling, atol=1.2*self.noise/quantize)

    def checkCompressedMaskedImage(self, image, imageOptions, maskOptions, varianceOptions, atol=0.0):
        """Check that compression works on a MaskedImage

        Parameters
        ----------
        image : `lsst.afw.image.MaskedImage`
            MaskedImage to compress.
        imageOptions, maskOptions, varianceOptions : `lsst.afw.fits.ImageWriteOptions`
            Parameters for writing (compression and scaling) the image, mask
            and variance planes.
        atol : `float`
            Absolute tolerance for comparing unpersisted image.
        """
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
        """Check that we can compress a MaskedImage and Exposure

        Parameters
        ----------
        imageOptions, maskOptions, varianceOptions : `lsst.afw.fits.ImageWriteOptions`
            Parameters for writing (compression and scaling) the image, mask
            and variance planes.
        atol : `float`
            Absolute tolerance for comparing unpersisted image.
        """
        image = lsst.afw.image.makeMaskedImage(self.makeImage(lsst.afw.image.ImageF),
                                               self.makeMask(), self.makeImage(lsst.afw.image.ImageF))
        self.checkCompressedMaskedImage(image, imageOptions, maskOptions, varianceOptions, atol=atol)
        exp = lsst.afw.image.makeExposure(image)
        self.checkCompressedMaskedImage(exp, imageOptions, maskOptions, varianceOptions, atol=atol)

    def testMaskedImage(self):
        """Test compression of MaskedImage

        We test lossless, lossy cfitsio and lossy LSST compression.
        """
        # Lossless
        lossless = lsst.afw.fits.ImageCompressionOptions(ImageCompressionOptions.GZIP_SHUFFLE)
        options = lsst.afw.fits.ImageWriteOptions(lossless)
        self.checkMaskedImage(options, options, options, atol=0.0)

        # Lossy cfitsio compression
        quantize = 4.0
        cfitsio = lsst.afw.fits.ImageCompressionOptions(ImageCompressionOptions.GZIP_SHUFFLE, True, quantize)
        imageOptions = lsst.afw.fits.ImageWriteOptions(cfitsio)
        maskOptions = lsst.afw.fits.ImageWriteOptions(lossless)
        self.checkMaskedImage(imageOptions, maskOptions, imageOptions, atol=1.25*self.noise/quantize)

        # Lossy our compression
        quantize = 4.0
        compression = lsst.afw.fits.ImageCompressionOptions(ImageCompressionOptions.GZIP_SHUFFLE, True, 0.0)
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
