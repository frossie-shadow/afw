#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

"""
Tests for Images

Run with:
   ./Image.py
or
   python
   >>> import Image; Image.run()
"""

from __future__ import absolute_import, division, print_function
import os.path
import shutil
import tempfile
import unittest

from builtins import range
import numpy as np

import lsst.utils
import lsst.utils.tests
import lsst.pex.exceptions
import lsst.daf.base
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.display.ds9 as ds9
import lsst.pex.exceptions as pexExcept

try:
    afwdataDir = lsst.utils.getPackageDir("afwdata")
except pexExcept.NotFoundError:
    afwdataDir = None

try:
    type(display)
except NameError:
    display = False


def makeRampImage(width, height, imgClass=afwImage.ImageF):
    """Make a ramp image of the specified size and image class

    Values start from 0 at the lower left corner and increase by 1 along rows
    """
    im = imgClass(width, height)
    val = 0
    for yInd in range(height):
        for xInd in range(width):
            im.set(xInd, yInd, val)
            val += 1
    return im


class ImageTestCase(lsst.utils.tests.TestCase):
    """A test case for Image"""

    def setUp(self):
        np.random.seed(1)
        self.val1, self.val2 = 10, 100
        self.image1 = afwImage.ImageF(afwGeom.ExtentI(100, 200))
        self.image1.set(self.val1)
        self.image2 = afwImage.ImageF(self.image1.getDimensions())
        self.image2.set(self.val2)
        self.function = afwMath.PolynomialFunction2D(2)
        self.function.setParameters(
            list(range(self.function.getNParameters())))

    def tearDown(self):
        del self.image1
        del self.image2
        del self.function

    def testArrays(self):
        for cls in (afwImage.ImageU, afwImage.ImageI, afwImage.ImageF, afwImage.ImageD):
            image1 = cls(afwGeom.Extent2I(5, 6))
            array1 = image1.getArray()
            self.assertEqual(array1.shape[0], image1.getHeight())
            self.assertEqual(array1.shape[1], image1.getWidth())
            image2 = cls(array1, False)
            self.assertEqual(array1.shape[0], image2.getHeight())
            self.assertEqual(array1.shape[1], image2.getWidth())
            image3 = afwImage.makeImageFromArray(array1)
            self.assertEqual(array1.shape[0], image2.getHeight())
            self.assertEqual(array1.shape[1], image2.getWidth())
            self.assertEqual(type(image3), cls)
            array2 = image1.array
            np.testing.assert_array_equal(array1, array2)
            array1[:, :] = np.random.uniform(low=0, high=10, size=array1.shape)
            for j in range(image1.getHeight()):
                for i in range(image1.getWidth()):
                    self.assertEqual(image1.get(i, j), array1[j, i])
                    self.assertEqual(image2.get(i, j), array1[j, i])
            array3 = np.random.uniform(low=0, high=10,
                                       size=array1.shape).astype(array1.dtype)
            image1.array[:] = array3
            np.testing.assert_array_equal(array1, array3)
            image1.array[2:4, 3:] = 10
            np.testing.assert_array_equal(array1[2:4, 3:], 10)
            array4 = image1.array.copy()
            array4 += 5
            image1.array += 5
            np.testing.assert_array_equal(image1.array, array4)

    def testInitializeImages(self):
        val = 666
        for ctor in (afwImage.ImageU, afwImage.ImageI, afwImage.ImageF, afwImage.ImageD):
            im = ctor(10, 10, val)
            self.assertEqual(im.get(0, 0), val)

            im2 = ctor(afwGeom.Box2I(afwGeom.Point2I(0, 0),
                                     afwGeom.Extent2I(10, 10)), val)
            self.assertEqual(im2.get(0, 0), val)

    def testSetGetImages(self):
        self.assertEqual(self.image1.get(0, 0), self.val1)

    def testGetSet0Images(self):
        self.assertEqual(self.image1.get0(0, 0), self.val1)
        self.image1.setXY0(3, 4)
        self.assertEqual(self.image1.get0(3, 4), self.val1)

        def f1():
            return self.image1.get0(0, 0)
        self.assertRaises(lsst.pex.exceptions.LengthError, f1)
        self.image1.set(0, 0, 42.)
        self.assertEqual(self.image1.get0(3, 4), 42.)
        self.image1.set0(3, 4, self.val1)
        self.assertEqual(self.image1.get0(3, 4), self.val1)
        self.assertEqual(self.image1.get(0, 0), self.val1)

    def testAllocateLargeImages(self):
        """Try to allocate a Very large image"""
        bbox = afwGeom.BoxI(afwGeom.PointI(-1 << 30, -1 << 30),
                            afwGeom.PointI(1 << 30, 1 << 30))

        def tst():
            afwImage.ImageF(bbox)

        self.assertRaises(lsst.pex.exceptions.LengthError, tst)

    def testAddImages(self):
        self.image2 += self.image1
        self.image1 += self.val1

        self.assertEqual(self.image1.get(0, 0), 2*self.val1)
        self.assertEqual(self.image2.get(0, 0), self.val1 + self.val2)

        self.image1.set(self.val1)
        self.image1 += self.function

        for j in range(self.image1.getHeight()):
            for i in range(self.image1.getWidth()):
                self.assertEqual(self.image1.get(i, j),
                                 self.val1 + self.function(i, j))

    def testAssignWithBBox(self):
        """Test assign(rhs, bbox) with non-empty bbox
        """
        for xy0 in (afwGeom.Point2I(*val) for val in (
            (0, 0),
            (-100, 120),  # an arbitrary value that is off the image
        )):
            destImDim = afwGeom.Extent2I(5, 4)
            srcImDim = afwGeom.Extent2I(3, 2)
            destIm = afwImage.ImageF(destImDim)
            destIm.setXY0(xy0)
            srcIm = makeRampImage(*srcImDim)
            srcIm.setXY0(55, -33)  # an arbitrary value that should be ignored
            self.assertRaises(Exception, destIm.set, srcIm)  # size mismatch

            for validMin in (afwGeom.Point2I(*val) for val in (
                (0, 0),
                (2, 0),
                (0, 1),
                (1, 2),
            )):
                # None to omit the argument
                for origin in (None, afwImage.PARENT, afwImage.LOCAL):
                    destIm[:] = -1.0
                    bbox = afwGeom.Box2I(validMin, srcIm.getDimensions())
                    if origin != afwImage.LOCAL:
                        bbox.shift(afwGeom.Extent2I(xy0))
                    if origin is None:
                        destIm.assign(srcIm, bbox)
                        destImView = afwImage.ImageF(destIm, bbox)
                    else:
                        destIm.assign(srcIm, bbox, origin)
                        destImView = afwImage.ImageF(destIm, bbox, origin)
                    self.assertFloatsEqual(
                        destImView.getArray(), srcIm.getArray())
                    numPixNotAssigned = (destImDim[0] * destImDim[1]) - \
                        (srcImDim[0] * srcImDim[1])
                    self.assertEqual(
                        np.sum(destIm.getArray() < -0.5), numPixNotAssigned)

            for badMin in (afwGeom.Point2I(*val) + afwGeom.Extent2I(xy0) for val in (
                (-1, 0),
                (3, 0),
                (0, -1),
                (1, 3),
            )):
                # None to omit the argument
                for origin in (None, afwImage.PARENT, afwImage.LOCAL):
                    bbox = afwGeom.Box2I(badMin, srcIm.getDimensions())
                    if origin != afwImage.LOCAL:
                        bbox.shift(afwGeom.Extent2I(xy0))
                    if origin is None:
                        self.assertRaises(Exception, destIm.set, srcIm, bbox)
                    else:
                        self.assertRaises(
                            Exception, destIm.set, srcIm, bbox, origin)

    def testAssignWithoutBBox(self):
        """Test assign(rhs, [bbox]) with an empty bbox and with no bbox specified; both set all pixels
        """
        for xy0 in (afwGeom.Point2I(*val) for val in (
            (0, 0),
            (-100, 120),  # an arbitrary value that is off the image
        )):
            destImDim = afwGeom.Extent2I(5, 4)
            destIm = afwImage.ImageF(destImDim)
            destIm.setXY0(xy0)
            srcIm = makeRampImage(*destImDim)
            srcIm.setXY0(55, -33)  # an arbitrary value that should be ignored

            destIm[:] = -1.0
            destIm.assign(srcIm)
            self.assertFloatsEqual(destIm.getArray(), srcIm.getArray())

            destIm[:] = -1.0
            destIm.assign(srcIm, afwGeom.Box2I())
            self.assertFloatsEqual(destIm.getArray(), srcIm.getArray())

    def testBoundsChecking(self):
        """Check that pixel indexes are checked in python"""
        tsts = []

        def tst():
            self.image1.get(-1, 0)
        tsts.append(tst)

        def tst():
            self.image1.get(0, -1)
        tsts.append(tst)

        def tst():
            self.image1.get(self.image1.getWidth(), 0)
        tsts.append(tst)

        def tst():
            self.image1.get(0, self.image1.getHeight())
        tsts.append(tst)

        for tst in tsts:
            self.assertRaises(lsst.pex.exceptions.LengthError, tst)

    def testAddScaledImages(self):
        c = 10.0
        self.image1.scaledPlus(c, self.image2)

        self.assertEqual(self.image1.get(0, 0), self.val1 + c*self.val2)

    def testSubtractImages(self):
        self.image2 -= self.image1
        self.image1 -= self.val1

        self.assertEqual(self.image1.get(0, 0), 0)
        self.assertEqual(self.image2.get(0, 0), self.val2 - self.val1)

        self.image1.set(self.val1)
        self.image1 -= self.function

        for j in range(self.image1.getHeight()):
            for i in range(self.image1.getWidth()):
                self.assertEqual(self.image1.get(i, j),
                                 self.val1 - self.function(i, j))

    def testArithmeticImagesMismatch(self):
        "Test arithmetic operations on Images of different sizes"
        i1 = afwImage.ImageF(100, 100)
        i1.set(100)
        i2 = afwImage.ImageF(10, 10)
        i2.set(10)

        def tst1(i1, i2):
            i1 -= i2

        def tst2(i1, i2):
            i1.scaledMinus(1.0, i2)

        def tst3(i1, i2):
            i1 += i2

        def tst4(i1, i2):
            i1.scaledPlus(1.0, i2)

        def tst5(i1, i2):
            i1 *= i2

        def tst6(i1, i2):
            i1.scaledMultiplies(1.0, i2)

        def tst7(i1, i2):
            i1 /= i2

        def tst8(i1, i2):
            i1.scaledDivides(1.0, i2)

        tsts12 = [tst1, tst3, tst5, tst7]
        for tst in tsts12:
            self.assertRaises(lsst.pex.exceptions.LengthError, tst, i1, i2)

        tsts21 = [tst2, tst4, tst6, tst8]
        for tst in tsts21:
            self.assertRaises(lsst.pex.exceptions.LengthError, tst, i2, i1)

    def testSubtractScaledImages(self):
        c = 10.0
        self.image1.scaledMinus(c, self.image2)

        self.assertEqual(self.image1.get(0, 0), self.val1 - c*self.val2)

    def testMultiplyImages(self):
        self.image2 *= self.image1
        self.image1 *= self.val1

        self.assertEqual(self.image1.get(0, 0), self.val1*self.val1)
        self.assertEqual(self.image2.get(0, 0), self.val2*self.val1)

    def testMultiplesScaledImages(self):
        c = 10.0
        self.image1.scaledMultiplies(c, self.image2)

        self.assertEqual(self.image1.get(0, 0), self.val1 * c*self.val2)

    def testDivideImages(self):
        self.image2 /= self.image1
        self.image1 /= self.val1

        self.assertEqual(self.image1.get(0, 0), 1)
        self.assertEqual(self.image2.get(0, 0), self.val2/self.val1)

    def testDividesScaledImages(self):
        c = 10.0
        self.image1.scaledDivides(c, self.image2)

        self.assertAlmostEqual(self.image1.get(0, 0), self.val1/(c*self.val2))

    def testCopyConstructors(self):
        dimage = afwImage.ImageF(self.image1, True)  # deep copy
        simage = afwImage.ImageF(self.image1)  # shallow copy

        self.image1 += 2                # should only change dimage
        self.assertEqual(dimage.get(0, 0), self.val1)
        self.assertEqual(simage.get(0, 0), self.val1 + 2)

    def testGeneralisedCopyConstructors(self):
        # these are generalised (templated) copy constructors in C++
        imageU = self.image1.convertU()
        imageF = imageU.convertF()
        imageD = imageF.convertD()

        self.assertEqual(imageU.get(0, 0), self.val1)
        self.assertEqual(imageF.get(0, 0), self.val1)
        self.assertEqual(imageD.get(0, 0), self.val1)

    def checkImgPatch(self, img, x0=0, y0=0):
        """Check that a patch of an image is correct; origin of patch is at (x0, y0)"""

        self.assertEqual(img.get(x0 - 1, y0 - 1), self.val1)
        self.assertEqual(img.get(x0, y0), 666)
        self.assertEqual(img.get(x0 + 3, y0), self.val1)
        self.assertEqual(img.get(x0, y0 + 1), 666)
        self.assertEqual(img.get(x0 + 3, y0 + 1), self.val1)
        self.assertEqual(img.get(x0, y0 + 2), self.val1)

    def testOrigin(self):
        """Check that we can set and read the origin"""

        im = afwImage.ImageF(10, 20)
        x0 = y0 = 0

        self.assertEqual(im.getX0(), x0)
        self.assertEqual(im.getY0(), y0)
        self.assertEqual(im.getXY0(), afwGeom.Point2I(x0, y0))

        x0, y0 = 3, 5
        im.setXY0(x0, y0)
        self.assertEqual(im.getX0(), x0)
        self.assertEqual(im.getY0(), y0)
        self.assertEqual(im.getXY0(), afwGeom.Point2I(x0, y0))

        x0, y0 = 30, 50
        im.setXY0(afwGeom.Point2I(x0, y0))
        self.assertEqual(im.getX0(), x0)
        self.assertEqual(im.getY0(), y0)
        self.assertEqual(im.getXY0(), afwGeom.Point2I(x0, y0))

    def testSubimages(self):
        simage1 = afwImage.ImageF(
            self.image1,
            afwGeom.Box2I(afwGeom.Point2I(1, 1), afwGeom.Extent2I(10, 5)),
            afwImage.LOCAL)

        simage = afwImage.ImageF(
            simage1,
            afwGeom.Box2I(afwGeom.Point2I(1, 1), afwGeom.Extent2I(3, 2)),
            afwImage.LOCAL
        )
        self.assertEqual(simage.getX0(), 2)
        self.assertEqual(simage.getY0(), 2)  # i.e. wrt self.image1

        image2 = afwImage.ImageF(simage.getDimensions())
        image2.set(666)
        simage[:] = image2
        del simage
        del image2

        self.checkImgPatch(self.image1, 2, 2)
        self.checkImgPatch(simage1, 1, 1)

    def testSubimages2(self):
        """Test subimages when we've played with the (x0, y0) value"""

        self.image1.set(9, 4, 888)

        simage1 = afwImage.ImageF(
            self.image1,
            afwGeom.Box2I(afwGeom.Point2I(1, 1), afwGeom.Extent2I(10, 5)),
            afwImage.LOCAL
        )
        # reset origin; doesn't affect pixel coordinate systems
        simage1.setXY0(afwGeom.Point2I(0, 0))

        simage = afwImage.ImageF(
            simage1,
            afwGeom.Box2I(afwGeom.Point2I(1, 1), afwGeom.Extent2I(3, 2)),
            afwImage.LOCAL
        )
        self.assertEqual(simage.getX0(), 1)
        self.assertEqual(simage.getY0(), 1)

        image2 = afwImage.ImageF(simage.getDimensions())
        image2.set(666)
        simage[:] = image2
        del simage
        del image2

        self.checkImgPatch(self.image1, 2, 2)
        self.checkImgPatch(simage1, 1, 1)

    def testBadSubimages(self):
        def tst():
            afwImage.ImageF(
                self.image1,
                afwGeom.Box2I(afwGeom.Point2I(1, -1), afwGeom.Extent2I(10, 5)),
                afwImage.LOCAL
            )

        self.assertRaises(lsst.pex.exceptions.LengthError, tst)

    def testImageInitialisation(self):
        dims = self.image1.getDimensions()
        factory = self.image1.Factory

        self.image1.set(666)

        del self.image1                 # tempt C++ to reuse the memory
        self.image1 = factory(dims)
        self.assertEqual(self.image1.get(10, 10), 0)

    def testImageSlices(self):
        """Test image slicing, which generate sub-images using Box2I under the covers"""
        im = afwImage.ImageF(10, 20)
        im[-1, :] = -5
        im[..., 18] = -5              # equivalent to im[:, 18]
        im[4, 10] = 10
        im[-3:, -2:] = 100
        im[-2, -2] = -10
        sim = im[1:4, 6:10]
        sim[:] = -1
        im[0:4, 0:4] = im[2:6, 8:12]

        if display:
            ds9.mtv(im)

        self.assertEqual(im.get(0, 6), 0)
        self.assertEqual(im.get(9, 15), -5)
        self.assertEqual(im.get(5, 18), -5)
        self.assertEqual(im.get(6, 17), 0)
        self.assertEqual(im.get(7, 18), 100)
        self.assertEqual(im.get(9, 19), 100)
        self.assertEqual(im.get(8, 18), -10)
        self.assertEqual(im.get(1, 6), -1)
        self.assertEqual(im.get(3, 9), -1)
        self.assertEqual(im.get(4, 10), 10)
        self.assertEqual(im.get(4, 9), 0)
        self.assertEqual(im.get(2, 2), 10)
        self.assertEqual(im.get(0, 0), -1)

    def testImageSlicesOrigin(self):
        """Test image slicing, which generate sub-images using Box2I under the covers"""
        im = afwImage.ImageF(10, 20)
        im.setXY0(50, 100)
        im[-1, :, afwImage.PARENT] = -5
        im[..., 118, afwImage.PARENT] = -5              # equivalent to im[:, 118]
        im[54, 110, afwImage.PARENT] = 10
        im[-3:, -2:, afwImage.PARENT] = 100
        im[-2, -2, afwImage.PARENT] = -10
        sim = im[51:54, 106:110, afwImage.PARENT]
        sim[:] = -1
        im[50:54, 100:104, afwImage.PARENT] = im[2:6, 8:12, afwImage.LOCAL]

        if display:
            ds9.mtv(im)

        self.assertEqual(im.get(0, 6), 0)
        self.assertEqual(im.get(9, 15), -5)
        self.assertEqual(im.get(5, 18), -5)
        self.assertEqual(im.get(6, 17), 0)
        self.assertEqual(im.get(7, 18), 100)
        self.assertEqual(im.get(9, 19), 100)
        self.assertEqual(im.get(8, 18), -10)
        self.assertEqual(im.get(1, 6), -1)
        self.assertEqual(im.get(3, 9), -1)
        self.assertEqual(im.get(4, 10), 10)
        self.assertEqual(im.get(4, 9), 0)
        self.assertEqual(im.get(2, 2), 10)
        self.assertEqual(im.get(0, 0), -1)

    def testImageSliceFromBox(self):
        """Test using a Box2I to index an Image"""
        im = afwImage.ImageF(10, 20)
        bbox = afwGeom.BoxI(afwGeom.PointI(1, 3), afwGeom.PointI(6, 9))
        im[bbox] = -1

        if display:
            ds9.mtv(im)

        self.assertEqual(im.get(0, 6), 0)
        self.assertEqual(im.get(1, 6), -1)
        self.assertEqual(im.get(3, 9), -1)

    def testImageSliceFromBoxOrigin(self):
        """Test using a Box2I to index an Image"""
        im = afwImage.ImageF(10, 20)
        im.setXY0(50, 100)
        bbox = afwGeom.BoxI(afwGeom.PointI(51, 103), afwGeom.ExtentI(6, 7))
        im[bbox, afwImage.PARENT] = -1

        if display:
            ds9.mtv(im)

        self.assertEqual(im.get(0, 6), 0)
        self.assertEqual(im.get(1, 6), -1)
        self.assertEqual(im.get(3, 9), -1)

    def testConversionToScalar(self):
        """Test that 1-pixel images can be converted to scalars"""
        self.assertEqual(int(afwImage.ImageI(1, 1)), 0.0)
        self.assertEqual(float(afwImage.ImageI(1, 1)), 0.0)

        im = afwImage.ImageF(10, 20)
        im.set(666)

        self.assertEqual(float(im[0, 0]), 666)
        self.assertEqual(int(im[0, 0]), 666)

        # only single pixel images may be converted
        self.assertRaises(TypeError, int, im)
        # only single pixel images may be converted
        self.assertRaises(TypeError, float, im)

    def testClone(self):
        """Test that clone works properly"""
        im = afwImage.ImageF(10, 20)
        im[0, 0] = 100

        im2 = im.clone()                # check that clone with no arguments makes a deep copy
        self.assertEqual(im.getDimensions(), im2.getDimensions())
        self.assertEqual(im.get(0, 0), im2.get(0, 0))
        im2[0, 0] += 100
        self.assertNotEqual(im.get(0, 0), im2.get(0, 0))  # so it's a deep copy

        im2 = im[0:3, 0:5].clone()  # check that we can slice-then-clone
        self.assertEqual(im2.getDimensions(), afwGeom.ExtentI(3, 5))
        self.assertEqual(im.get(0, 0), im2.get(0, 0))
        im2[0, 0] += 10
        # equivalent to im.get(0, 0) etc.
        self.assertNotEqual(float(im[0, 0]), float(im2[0, 0]))


class DecoratedImageTestCase(lsst.utils.tests.TestCase):
    """A test case for DecoratedImage"""

    def setUp(self):
        np.random.seed(1)
        self.val1, self.val2 = 10, 100
        self.width, self.height = 200, 100
        self.dimage1 = afwImage.DecoratedImageF(
            afwGeom.Extent2I(self.width, self.height)
        )
        self.dimage1.getImage().set(self.val1)

        if afwdataDir is not None:
            self.fileForMetadata = os.path.join(
                afwdataDir, "data", "small_MI.fits")
            self.trueMetadata = {"RELHUMID": 10.69}

    def tearDown(self):
        del self.dimage1

    def testCreateDecoratedImage(self):
        self.assertEqual(self.dimage1.getWidth(), self.width)
        self.assertEqual(self.dimage1.getHeight(), self.height)
        self.assertEqual(self.dimage1.getImage().get(0, 0), self.val1)

    def testCreateDecoratedImageFromImage(self):
        image = afwImage.ImageF(afwGeom.Extent2I(self.width, self.height))
        image[:] = self.dimage1.getImage()

        dimage = afwImage.DecoratedImageF(image)
        self.assertEqual(dimage.getWidth(), self.width)
        self.assertEqual(dimage.getHeight(), self.height)
        self.assertEqual(dimage.getImage().get(0, 0), self.val1)

    def testCopyConstructors(self):
        dimage = afwImage.DecoratedImageF(self.dimage1, True)  # deep copy
        self.dimage1.getImage().set(0, 0, 1 + 2*self.val1)
        self.assertEqual(dimage.getImage().get(0, 0), self.val1)

        dimage = afwImage.DecoratedImageF(self.dimage1)  # shallow copy
        self.dimage1.getImage().set(0, 0, 1 + 2*self.val1)
        self.assertNotEqual(dimage.getImage().get(0, 0), self.val1)

    @unittest.skipIf(afwdataDir is None, "afwdata not setup")
    def testReadFits(self):
        """Test reading FITS files"""

        hdus = {}
        hdus["img"] = 1  # an S16 fits HDU
        hdus["msk"] = 2  # an U8 fits HDU
        hdus["var"] = 3  # an F32 fits HDU

        # read as unsigned short
        imgU = afwImage.DecoratedImageU(self.fileForMetadata, hdus["img"])
        # read as float
        imgF = afwImage.DecoratedImageF(self.fileForMetadata, hdus["img"])

        self.assertEqual(imgU.getHeight(), 256)
        self.assertEqual(imgF.getImage().getWidth(), 256)
        self.assertEqual(imgU.getImage().get(0, 0), imgF.getImage().get(0, 0))
        #
        # Check the metadata
        #
        meta = self.trueMetadata
        for k in meta.keys():
            self.assertEqual(imgU.getMetadata().getAsDouble(k), meta[k])
            self.assertEqual(imgF.getMetadata().getAsDouble(k), meta[k])
        #
        # Read an F32 image
        #
        # read as unsigned short
        varU = afwImage.DecoratedImageF(self.fileForMetadata, hdus["var"])
        # read as float
        varF = afwImage.DecoratedImageF(self.fileForMetadata, hdus["var"])

        self.assertEqual(varU.getHeight(), 256)
        self.assertEqual(varF.getImage().getWidth(), 256)
        self.assertEqual(varU.getImage().get(0, 0), varF.getImage().get(0, 0))
        #
        # Read a char image
        #
        maskImg = afwImage.DecoratedImageU(
            self.fileForMetadata, hdus["msk"]).getImage()  # read a char file

        self.assertEqual(maskImg.getHeight(), 256)
        self.assertEqual(maskImg.getWidth(), 256)
        self.assertEqual(maskImg.get(0, 0), 1)
        #
        # Read a U16 image
        #
        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            imgU.writeFits(tmpFile)

            afwImage.DecoratedImageF(tmpFile)  # read as unsigned short

    @unittest.skipIf(afwdataDir is None, "afwdata not setup")
    def testWriteFits(self):
        """Test writing FITS files"""

        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            if self.fileForMetadata:
                imgU = afwImage.DecoratedImageF(self.fileForMetadata)
            else:
                imgU = afwImage.DecoratedImageF()

            self.dimage1.writeFits(tmpFile, imgU.getMetadata())
            #
            # Read it back
            #
            rimage = afwImage.DecoratedImageF(tmpFile)

            self.assertEqual(self.dimage1.getImage().get(0, 0),
                             rimage.getImage().get(0, 0))
            #
            # Check that we wrote (and read) the metadata successfully
            if self.fileForMetadata:
                meta = self.trueMetadata
                for k in meta.keys():
                    self.assertEqual(
                        rimage.getMetadata().getAsDouble(k), meta[k])

    def testReadWriteXY0(self):
        """Test that we read and write (X0, Y0) correctly"""
        im = afwImage.ImageF(afwGeom.Extent2I(10, 20))

        x0, y0 = 1, 2
        im.setXY0(x0, y0)
        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            im.writeFits(tmpFile)

            im2 = im.Factory(tmpFile)

            self.assertEqual(im2.getX0(), x0)
            self.assertEqual(im2.getY0(), y0)

    @unittest.skipIf(afwdataDir is None, "afwdata not setup")
    def testReadMetadata(self):
        im = afwImage.DecoratedImageF(self.fileForMetadata)

        meta = afwImage.readMetadata(self.fileForMetadata)
        self.assertIn("NAXIS1", meta.names())
        self.assertEqual(im.getWidth(), meta.get("NAXIS1"))
        self.assertEqual(im.getHeight(), meta.get("NAXIS2"))

    def testTicket1040(self):
        """ How to repeat from #1040"""
        image = afwImage.ImageD(afwGeom.Extent2I(6, 6))
        image.set(2, 2, 100)

        bbox = afwGeom.Box2I(afwGeom.Point2I(1, 1), afwGeom.Extent2I(5, 5))
        subImage = image.Factory(image, bbox)
        subImageF = subImage.convertFloat()

        if display:
            ds9.mtv(subImage, frame=0, title="subImage")
            ds9.mtv(subImageF, frame=1, title="converted subImage")

        self.assertEqual(subImage.get(1, 1), subImageF.get(1, 1))

    def testDM882(self):
        """Test that we can write a dotted header unit to a FITS file. See DM-882."""
        self.dimage1.getMetadata().add("A.B.C.D", 12345)
        tempdir = tempfile.mkdtemp()
        testfile = os.path.join(tempdir, "test.fits")
        try:
            self.dimage1.writeFits(testfile)
            meta = afwImage.readMetadata(testfile)
            self.assertEqual(meta.get("A.B.C.D"), 12345)
        finally:
            shutil.rmtree(tempdir)

    def testLargeImage(self):
        """Test that creating an extremely large image raises, rather than segfaulting. DM-89, -527."""
        for imtype in (afwImage.ImageD, afwImage.ImageF, afwImage.ImageI, afwImage.ImageU):
            self.assertRaises(lsst.pex.exceptions.LengthError,
                              imtype, 60000, 60000)


def printImg(img):
    print("%4s " % "", end=' ')
    for c in range(img.getWidth()):
        print("%7d" % c, end=' ')
    print()

    for r in range(img.getHeight() - 1, -1, -1):
        print("%4d " % r, end=' ')
        for c in range(img.getWidth()):
            print("%7.1f" % float(img.get(c, r)), end=' ')
        print()


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
