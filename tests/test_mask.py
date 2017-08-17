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
Tests for Masks

Run with:
   python Mask.py
or
   python
   >>> import Mask; Mask.run()
"""
from __future__ import absolute_import, division, print_function
import os.path
import unittest

from builtins import range
import numpy as np

import lsst.utils
import lsst.utils.tests as utilsTests
import lsst.pex.exceptions as pexExcept
import lsst.daf.base
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.display.ds9 as ds9


try:
    type(display)
except NameError:
    display = False

try:
    afwdataDir = lsst.utils.getPackageDir("afwdata")
except pexExcept.NotFoundError:
    afwdataDir = None


def showMaskDict(d=None, msg=None):
    if not d:
        d = afwImage.Mask(0, 0)
        if not msg:
            msg = "default"

    try:
        d = d.getMaskPlaneDict()
    except AttributeError:
        pass

    if msg:
        print("%-15s" % msg, end=' ')
    print(sorted([(d[p], p) for p in d]))


class MaskTestCase(utilsTests.TestCase):
    """A test case for Mask"""

    def setUp(self):
        np.random.seed(1)
        self.Mask = afwImage.Mask[afwImage.MaskPixel]

        # Store the default mask planes for later use
        maskPlaneDict = self.Mask().getMaskPlaneDict()
        self.defaultMaskPlanes = sorted(
            maskPlaneDict, key=maskPlaneDict.__getitem__)

        # reset so tests will be deterministic
        self.Mask.clearMaskPlaneDict()
        for p in ("BAD", "SAT", "INTRP", "CR", "EDGE"):
            self.Mask.addMaskPlane(p)

        self.BAD = self.Mask.getPlaneBitMask("BAD")
        self.CR = self.Mask.getPlaneBitMask("CR")
        self.EDGE = self.Mask.getPlaneBitMask("EDGE")

        self.val1 = self.BAD | self.CR
        self.val2 = self.val1 | self.EDGE

        self.mask1 = afwImage.Mask(100, 200)
        self.mask1.set(self.val1)
        self.mask2 = afwImage.Mask(self.mask1.getDimensions())
        self.mask2.set(self.val2)

        # TBD: #DM-609 this should be refactored to use @unittest.skipif checks
        # for afwData above the tests that need it.
        if afwdataDir is not None:
            self.maskFile = os.path.join(afwdataDir, "data", "small_MI.fits")
            # Below: what to expect from the mask plane in the above data.
            # For some tests, it is left-shifted by some number of mask planes.
            self.expect = np.zeros((256, 256), dtype='i8')
            self.expect[:, 0] = 1

    def tearDown(self):
        del self.mask1
        del self.mask2
        # Reset the mask plane to the default
        self.Mask.clearMaskPlaneDict()
        for p in self.defaultMaskPlanes:
            self.Mask.addMaskPlane(p)

    def testArrays(self):
        # could use Mask(5, 6) but check extent(5, 6) form too
        image1 = afwImage.Mask(afwGeom.ExtentI(5, 6))
        array1 = image1.getArray()
        self.assertEqual(array1.shape[0], image1.getHeight())
        self.assertEqual(array1.shape[1], image1.getWidth())
        image2 = afwImage.Mask(array1, False)
        self.assertEqual(array1.shape[0], image2.getHeight())
        self.assertEqual(array1.shape[1], image2.getWidth())
        image3 = afwImage.makeMaskFromArray(array1)
        self.assertEqual(array1.shape[0], image2.getHeight())
        self.assertEqual(array1.shape[1], image2.getWidth())
        self.assertEqual(type(image3), afwImage.Mask[afwImage.MaskPixel])
        array1[:, :] = np.random.uniform(low=0, high=10, size=array1.shape)
        self.assertMasksEqual(image1, array1)
        array2 = image1.array
        np.testing.assert_array_equal(image1.array, array2)
        array3 = np.random.uniform(low=0, high=10,
                                   size=array1.shape).astype(array1.dtype)
        image1.array = array3
        np.testing.assert_array_equal(array1, array3)

    def testInitializeMasks(self):
        val = 0x1234
        msk = afwImage.Mask(afwGeom.ExtentI(10, 10), val)
        self.assertEqual(msk.get(0, 0), val)

    def testSetGetMasks(self):
        self.assertEqual(self.mask1.get(0, 0), self.val1)

    def testOrMasks(self):
        self.mask2 |= self.mask1
        self.mask1 |= self.val2

        self.assertMasksEqual(self.mask1, self.mask2)
        expect = np.empty_like(self.mask1.getArray())
        expect[:] = self.val2 | self.val1
        self.assertMasksEqual(self.mask1, expect)

    def testAndMasks(self):
        self.mask2 &= self.mask1
        self.mask1 &= self.val2

        self.assertMasksEqual(self.mask1, self.mask2)
        expect = np.empty_like(self.mask1.getArray())
        expect[:] = self.val1 & self.val2
        self.assertMasksEqual(self.mask1, expect)

    def testXorMasks(self):
        self.mask2 ^= self.mask1
        self.mask1 ^= self.val2

        self.assertMasksEqual(self.mask1, self.mask2)
        expect = np.empty_like(self.mask1.getArray())
        expect[:] = self.val1 ^ self.val2
        self.assertMasksEqual(self.mask1, expect)

    def testLogicalMasksMismatch(self):
        "Test logical operations on Masks of different sizes"
        i1 = afwImage.Mask(afwGeom.ExtentI(100, 100))
        i1.set(100)
        i2 = afwImage.Mask(afwGeom.ExtentI(10, 10))
        i2.set(10)

        with self.assertRaises(lsst.pex.exceptions.LengthError):
            i1 |= i2

        with self.assertRaises(lsst.pex.exceptions.LengthError):
            i1 &= i2

        with self.assertRaises(lsst.pex.exceptions.LengthError):
            i1 ^= i2

    def testMaskPlanes(self):
        planes = self.Mask().getMaskPlaneDict()
        self.assertEqual(len(planes), self.Mask.getNumPlanesUsed())

        for k in sorted(planes.keys()):
            self.assertEqual(planes[k], self.Mask.getMaskPlane(k))

    def testCopyConstructors(self):
        dmask = afwImage.Mask(self.mask1, True)  # deep copy
        smask = afwImage.Mask(self.mask1)  # shallow copy

        self.mask1 |= 32767             # should only change smask
        temp = np.zeros_like(self.mask1.getArray()) | self.val1
        self.assertMasksEqual(dmask, temp)
        self.assertMasksEqual(smask, self.mask1)
        self.assertMasksEqual(smask, temp | 32767)

    def testSubmasks(self):
        smask = afwImage.Mask(self.mask1,
                              afwGeom.Box2I(afwGeom.Point2I(1, 1),
                                            afwGeom.ExtentI(3, 2)),
                              afwImage.LOCAL)
        mask2 = afwImage.Mask(smask.getDimensions())

        mask2.set(666)
        smask[:] = mask2

        del smask
        del mask2

        self.assertEqual(self.mask1.get(0, 0), self.val1)
        self.assertEqual(self.mask1.get(1, 1), 666)
        self.assertEqual(self.mask1.get(4, 1), self.val1)
        self.assertEqual(self.mask1.get(1, 2), 666)
        self.assertEqual(self.mask1.get(4, 2), self.val1)
        self.assertEqual(self.mask1.get(1, 3), self.val1)

    @unittest.skipIf(afwdataDir is None, "afwdata not setup")
    def testReadFits(self):
        nMaskPlanes0 = self.Mask.getNumPlanesUsed()
        # will shift any unrecognised mask planes into unused slots
        mask = self.Mask(self.maskFile, hdu=2)

        self.assertMasksEqual(mask, self.expect << nMaskPlanes0)

    @unittest.skipIf(afwdataDir is None, "afwdata not setup")
    def testReadFitsConform(self):
        hdu = 2
        mask = afwImage.Mask(self.maskFile, hdu, None,
                             afwGeom.Box2I(), afwImage.LOCAL, True)

        self.assertMasksEqual(mask, self.expect)

    @unittest.skipIf(afwdataDir is None, "afwdata not setup")
    def testWriteFits(self):
        nMaskPlanes0 = self.Mask.getNumPlanesUsed()
        mask = self.Mask(self.maskFile, hdu=2)

        self.assertMasksEqual(mask, self.expect << nMaskPlanes0)

        mask.clearMaskPlaneDict()

        with utilsTests.getTempFilePath(".fits") as tmpFile:
            mask.writeFits(tmpFile)

            # Read it back
            md = lsst.daf.base.PropertySet()
            rmask = self.Mask(tmpFile, 0, md)
            self.assertMasksEqual(mask, rmask)

            # Check that we wrote (and read) the metadata successfully
            mp_ = "MP_" if True else self.Mask.maskPlanePrefix()  # currently private
            for (k, v) in self.Mask().getMaskPlaneDict().items():
                self.assertEqual(md.get(mp_ + k), v)

    def testReadWriteXY0(self):
        """Test that we read and write (X0, Y0) correctly"""
        mask = afwImage.Mask(afwGeom.ExtentI(10, 20))

        x0, y0 = 1, 2
        mask.setXY0(x0, y0)
        with utilsTests.getTempFilePath(".fits") as tmpFile:
            mask.writeFits(tmpFile)

            mask2 = mask.Factory(tmpFile)

            self.assertEqual(mask2.getX0(), x0)
            self.assertEqual(mask2.getY0(), y0)

    def testMaskInitialisation(self):
        dims = self.mask1.getDimensions()
        factory = self.mask1.Factory

        self.mask1.set(666)

        del self.mask1                 # tempt C++ to reuse the memory
        self.mask1 = factory(dims)
        self.assertEqual(self.mask1.get(10, 10), 0)

        del self.mask1
        self.mask1 = factory(afwGeom.ExtentI(20, 20))
        self.assertEqual(self.mask1.get(10, 10), 0)

    def testBoundsChecking(self):
        """Check that pixel indexes are checked in python"""
        tsts = []

        def tst():
            self.mask1.get(-1, 0)
        tsts.append(tst)

        def tst():
            self.mask1.get(0, -1)
        tsts.append(tst)

        def tst():
            self.mask1.get(self.mask1.getWidth(), 0)
        tsts.append(tst)

        def tst():
            self.mask1.get(0, self.mask1.getHeight())
        tsts.append(tst)

        for tst in tsts:
            self.assertRaises(lsst.pex.exceptions.LengthError, tst)

    def testCtorWithPlaneDefs(self):
        """Test that we can create a Mask with a given MaskPlaneDict"""
        FOO, val = "FOO", 2
        mask = afwImage.Mask(100, 200, {FOO: val}
                              )
        mpd = mask.getMaskPlaneDict()
        self.assertIn(FOO, mpd.keys())
        self.assertEqual(mpd[FOO], val)

    def testImageSlices(self):
        """Test image slicing, which generate sub-images using Box2I under the covers"""
        im = afwImage.Mask(10, 20)
        im[-3:, -2:] = 0x4
        im[4, 10] = 0x2
        sim = im[1:4, 6:10]
        sim[:] = 0x8
        im[0:4, 0:4] = im[2:6, 8:12]

        if display:
            ds9.mtv(im)

        self.assertEqual(im.get(0, 6), 0)
        self.assertEqual(im.get(6, 17), 0)
        self.assertEqual(im.get(7, 18), 0x4)
        self.assertEqual(im.get(9, 19), 0x4)
        self.assertEqual(im.get(1, 6), 0x8)
        self.assertEqual(im.get(3, 9), 0x8)
        self.assertEqual(im.get(4, 10), 0x2)
        self.assertEqual(im.get(4, 9), 0)
        self.assertEqual(im.get(2, 2), 0x2)
        self.assertEqual(im.get(0, 0), 0x8)

    def testInterpret(self):
        """Interpretation of Mask values"""
        planes = self.Mask().getMaskPlaneDict()
        im = self.Mask(len(planes), 1)

        allBits = 0
        for i, p in enumerate(planes):
            bitmask = self.Mask.getPlaneBitMask(p)
            allBits |= bitmask
            self.assertEqual(im.interpret(bitmask), p)
            im.getArray()[0, i] = bitmask
            self.assertEqual(im.getAsString(i, 0), p)
        self.assertEqual(self.Mask.interpret(allBits),
                         ",".join(sorted(planes.keys())))


class OldMaskTestCase(unittest.TestCase):
    """A test case for Mask (based on Mask_1.cc); these are taken over from the DC2 fw tests
    and modified to run with the new (DC3) APIs"""

    def setUp(self):
        self.Mask = afwImage.Mask[afwImage.MaskPixel]           # the class

        self.testMask = self.Mask(afwGeom.Extent2I(300, 400), 0)

        # Store the default mask planes for later use
        maskPlaneDict = self.Mask().getMaskPlaneDict()
        self.defaultMaskPlanes = sorted(
            maskPlaneDict, key=maskPlaneDict.__getitem__)

        # reset so tests will be deterministic
        self.Mask.clearMaskPlaneDict()
        for p in ("CR", "BP"):
            self.Mask.addMaskPlane(p)

        self.region = afwGeom.Box2I(afwGeom.Point2I(100, 300),
                                    afwGeom.Extent2I(10, 40))
        self.subTestMask = self.Mask(
            self.testMask, self.region, afwImage.LOCAL)

        if False:
            self.pixelList = afwImage.listPixelCoord()
            for x in range(0, 300):
                for y in range(300, 400, 20):
                    self.pixelList.push_back(afwImage.PixelCoord(x, y))

    def tearDown(self):
        del self.testMask
        del self.subTestMask
        del self.region
        # Reset the mask plane to the default
        self.Mask.clearMaskPlaneDict()
        for p in self.defaultMaskPlanes:
            self.Mask.addMaskPlane(p)

    def testPlaneAddition(self):
        """Test mask plane addition"""

        nplane = self.testMask.getNumPlanesUsed()
        for p in ("XCR", "XBP"):
            self.assertEqual(self.Mask.addMaskPlane(p),
                             nplane, "Assigning plane %s" % (p))
            nplane += 1

        def pname(p):
            return "P%d" % p

        nextra = 8
        for p in range(0, nextra):
            self.Mask.addMaskPlane(pname(p))

        for p in range(0, nextra):
            self.testMask.removeAndClearMaskPlane(pname(p))

        self.assertEqual(
            nplane + nextra, self.Mask.getNumPlanesUsed(), "Adding and removing planes")

        for p in range(0, nextra):
            self.Mask.removeMaskPlane(pname(p))

        self.assertEqual(nplane, self.testMask.getNumPlanesUsed(),
                         "Adding and removing planes")

    def testMetadata(self):
        """Test mask plane metadata interchange with MaskPlaneDict"""
        #
        # Demonstrate that we can extract a MaskPlaneDict into metadata
        #
        metadata = lsst.daf.base.PropertySet()

        self.Mask.addMaskPlanesToMetadata(metadata)
        for (k, v) in self.Mask().getMaskPlaneDict().items():
            self.assertEqual(metadata.getInt("MP_%s" % k), v)
        #
        # Now add another plane to metadata and make it appear in the mask Dict, albeit
        # in general at another location (hence the getNumPlanesUsed call)
        #
        metadata.addInt("MP_" + "Whatever", self.Mask.getNumPlanesUsed())

        self.testMask.conformMaskPlanes(
            self.Mask.parseMaskPlaneMetadata(metadata))
        for (k, v) in self.Mask().getMaskPlaneDict().items():
            self.assertEqual(metadata.getInt("MP_%s" % k), v)

    def testPlaneOperations(self):
        """Test mask plane operations"""

        planes = self.Mask().getMaskPlaneDict()
        self.testMask.clearMaskPlane(planes['CR'])

        if False:
            for p in planes.keys():
                self.testMask.setMaskPlaneValues(planes[p], self.pixelList)

        # print "\nClearing mask"
        self.testMask.clearMaskPlane(planes['CR'])

    def testPlaneRemoval(self):
        """Test mask plane removal"""

        def checkPlaneBP():
            self.Mask.getMaskPlane("BP")

        testMask2 = self.Mask(self.testMask.getDimensions())
        self.testMask = self.Mask(self.testMask.getDimensions())
        self.testMask.removeAndClearMaskPlane("BP")

        testMask2.getMaskPlaneDict()

        # still present in default mask
        checkPlaneBP()
        # should still be in testMask2
        self.assertIn("BP", testMask2.getMaskPlaneDict())

        self.Mask.removeMaskPlane("BP")  # remove from default mask too

        self.assertRaises(pexExcept.InvalidParameterError, checkPlaneBP)

        self.assertRaises(pexExcept.InvalidParameterError,
                          lambda: self.Mask.removeMaskPlane("BP"))  # Plane is already removed
        self.assertRaises(pexExcept.InvalidParameterError,
                          lambda: self.testMask.removeMaskPlane("RHL gets names right"))
        #
        self.Mask.clearMaskPlaneDict()
        self.Mask.addMaskPlane("P0")
        self.Mask.addMaskPlane("P1")
        # a no-op -- readding a plane has no effect
        self.Mask.addMaskPlane("P1")
        #
        # Check that removing default mask planes doesn't affect pre-existing planes
        #
        msk = self.Mask()
        nmask = len(msk.getMaskPlaneDict())
        self.Mask.removeMaskPlane("P0")
        self.Mask.removeMaskPlane("P1")
        self.assertEqual(len(msk.getMaskPlaneDict()), nmask)
        del msk
        #
        # Check that removeAndClearMaskPlane can clear the default too
        #
        self.Mask.addMaskPlane("BP")
        self.testMask.removeAndClearMaskPlane("BP", True)

        self.assertRaises(pexExcept.InvalidParameterError, checkPlaneBP)

    def testInvalidPlaneOperations(self):
        """Test mask plane operations invalidated by Mask changes"""

        testMask3 = self.Mask(self.testMask.getDimensions())

        name = "Great Timothy"
        self.Mask.addMaskPlane(name)
        testMask3.removeAndClearMaskPlane(name)

        self.Mask.getMaskPlane(name)    # should be fine
        self.assertRaises(KeyError, lambda: testMask3.getMaskPlaneDict()[name])

        def tst():
            self.testMask |= testMask3

        self.assertRaises(pexExcept.RuntimeError, tst)

        # The dictionary should be back to the same state, so ...
        self.Mask.addMaskPlane(name)
        tst                             # ... assertion should not fail

        self.testMask.removeAndClearMaskPlane(name, True)
        self.Mask.addMaskPlane("Mario")  # takes name's slot
        self.Mask.addMaskPlane(name)

        self.assertRaises(pexExcept.RuntimeError, tst)

    def testInvalidPlaneOperations2(self):
        """Test mask plane operations invalidated by Mask changes"""

        testMask3 = self.Mask(self.testMask.getDimensions())

        name = "Great Timothy"
        name2 = "Our Boss"
        self.Mask.addMaskPlane(name)
        self.Mask.addMaskPlane(name2)
        # a description of the Mask's current dictionary
        oldDict = testMask3.getMaskPlaneDict()

        for n in (name, name2):
            self.testMask.removeAndClearMaskPlane(n, True)

        # added in opposite order to the planes in testMask3
        self.Mask.addMaskPlane(name2)
        self.Mask.addMaskPlane(name)

        self.assertNotEqual(self.testMask.getMaskPlaneDict()[
                            name], oldDict[name])

        def tst():
            self.testMask |= testMask3

        self.testMask.removeAndClearMaskPlane("BP")

        self.assertRaises(pexExcept.RuntimeError, tst)
        #
        # OK, that failed as it should.  Fixup the dictionaries and try again
        #
        self.Mask.addMaskPlane("BP")
        # convert testMask3 from oldDict to current default
        testMask3.conformMaskPlanes(oldDict)

        self.testMask |= testMask3      # shouldn't throw

    def testConformMaskPlanes(self):
        """Test conformMaskPlanes() when the two planes are actually the same"""

        testMask3 = self.Mask(self.testMask.getDimensions())

        name = "XXX"
        self.Mask.addMaskPlane(name)
        oldDict = testMask3.getMaskPlaneDict()
        # invalidates dictionary version
        testMask3.removeAndClearMaskPlane(name)

        testMask3.conformMaskPlanes(oldDict)

        self.testMask |= testMask3

    def testConformMaskPlanes2(self):
        """Test conformMaskPlanes() when the two planes are different"""

        testMask3 = afwImage.Mask(self.testMask.getDimensions())

        name1 = "Great Timothy"
        name2 = "Our Boss"
        p1 = self.Mask.addMaskPlane(name1)
        p2 = self.Mask.addMaskPlane(name2)
        oldDict = self.testMask.getMaskPlaneDict()

        testMask3.setMaskPlaneValues(p1, 0, 5, 0)
        testMask3.setMaskPlaneValues(p2, 0, 5, 1)

        if display:
            im = afwImage.ImageF(testMask3.getDimensions())
            # bug in ds9's Mask display; needs an Image first
            ds9.mtv(im)
            ds9.mtv(testMask3)

        self.assertEqual(testMask3.get(0, 0), testMask3.getPlaneBitMask(name1))
        self.assertEqual(testMask3.get(0, 1), testMask3.getPlaneBitMask(name2))

        self.testMask.removeAndClearMaskPlane(name1, True)
        self.testMask.removeAndClearMaskPlane(name2, True)
        self.Mask.addMaskPlane(name2)  # added in opposite order to testMask3
        self.Mask.addMaskPlane(name1)

        self.assertEqual(self.testMask.get(0, 0), 0)

        if display:
            ds9.mtv(im, frame=1)
            ds9.mtv(testMask3, frame=1)

        self.assertNotEqual(testMask3.get(0, 0),
                            testMask3.getPlaneBitMask(name1))
        self.assertNotEqual(testMask3.get(0, 1),
                            testMask3.getPlaneBitMask(name2))

        testMask3.conformMaskPlanes(oldDict)

        self.assertEqual(testMask3.get(0, 0), testMask3.getPlaneBitMask(name1))
        self.assertEqual(testMask3.get(0, 1), testMask3.getPlaneBitMask(name2))

        if display:
            ds9.mtv(im, frame=2)
            ds9.mtv(testMask3, frame=2)

        self.testMask |= testMask3


def printMaskPlane(mask, plane,
                   xrange=list(range(250, 300, 10)), yrange=list(range(300, 400, 20))):
    """Print parts of the specified plane of the mask"""

    xrange = list(range(min(xrange), max(xrange), 25))
    yrange = list(range(min(yrange), max(yrange), 25))

    for x in xrange:
        for y in yrange:
            if False:                   # mask(x,y) confuses swig
                print(x, y, mask(x, y), mask(x, y, plane))
            else:
                print(x, y, mask(x, y, plane))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
