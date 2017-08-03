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
Test lsst.afw.image.Exposure
"""

from __future__ import absolute_import, division, print_function
import os.path
import unittest

from builtins import range
import numpy as np

import lsst.utils
import lsst.utils.tests
import lsst.afw.image as afwImage
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.pex.exceptions as pexExcept
import lsst.pex.policy as pexPolicy
import lsst.afw.fits
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
from lsst.log import Log
from testTableArchivesLib import DummyPsf

Log.getLogger("afw.image.Mask").setLevel(Log.INFO)

try:
    dataDir = os.path.join(lsst.utils.getPackageDir("afwdata"), "data")
except pexExcept.NotFoundError:
    dataDir = None
else:
    InputMaskedImageName = "871034p_1_MI.fits"
    InputMaskedImageNameSmall = "small_MI.fits"
    InputImageNameSmall = "small"
    OutputMaskedImageName = "871034p_1_MInew.fits"

    currDir = os.path.abspath(os.path.dirname(__file__))
    inFilePath = os.path.join(dataDir, InputMaskedImageName)
    inFilePathSmall = os.path.join(dataDir, InputMaskedImageNameSmall)
    inFilePathSmallImage = os.path.join(dataDir, InputImageNameSmall)


@unittest.skipIf(dataDir is None, "afwdata not setup")
class ExposureTestCase(lsst.utils.tests.TestCase):
    """
    A test case for the Exposure Class
    """

    def setUp(self):
        maskedImage = afwImage.MaskedImageF(inFilePathSmall)
        maskedImageMD = afwImage.readMetadata(inFilePathSmall)

        self.smallExposure = afwImage.ExposureF(inFilePathSmall)
        self.width = maskedImage.getWidth()
        self.height = maskedImage.getHeight()
        self.wcs = afwGeom.SkyWcs(maskedImageMD)
        self.psf = DummyPsf(2.0)
        self.detector = DetectorWrapper().detector

        self.exposureBlank = afwImage.ExposureF()
        self.exposureMiOnly = afwImage.makeExposure(maskedImage)
        self.exposureMiWcs = afwImage.makeExposure(maskedImage, self.wcs)
        # n.b. the (100, 100, ...) form
        self.exposureCrWcs = afwImage.ExposureF(100, 100, self.wcs)
        # test with ExtentI(100, 100) too
        self.exposureCrOnly = afwImage.ExposureF(afwGeom.ExtentI(100, 100))

        afwImage.Filter.reset()
        afwImage.FilterProperty.reset()

        filterPolicy = pexPolicy.Policy()
        filterPolicy.add("lambdaEff", 470.0)
        afwImage.Filter.define(afwImage.FilterProperty("g", filterPolicy))

    def tearDown(self):
        del self.smallExposure
        del self.wcs
        del self.psf
        del self.detector

        del self.exposureBlank
        del self.exposureMiOnly
        del self.exposureMiWcs
        del self.exposureCrWcs
        del self.exposureCrOnly

    def testGetMaskedImage(self):
        """
        Test to ensure a MaskedImage can be obtained from each
        Exposure. An Exposure is required to have a MaskedImage,
        therefore each of the Exposures should return a MaskedImage.

        MaskedImage class should throw appropriate
        lsst::pex::exceptions::NotFound if the MaskedImage can not be
        obtained.
        """
        maskedImageBlank = self.exposureBlank.getMaskedImage()
        blankWidth = maskedImageBlank.getWidth()
        blankHeight = maskedImageBlank.getHeight()
        if blankWidth != blankHeight != 0:
            self.fail("%s = %s != 0" % (blankWidth, blankHeight))

        maskedImageMiOnly = self.exposureMiOnly.getMaskedImage()
        miOnlyWidth = maskedImageMiOnly.getWidth()
        miOnlyHeight = maskedImageMiOnly.getHeight()
        self.assertAlmostEqual(miOnlyWidth, self.width)
        self.assertAlmostEqual(miOnlyHeight, self.height)

        # NOTE: Unittests for Exposures created from a MaskedImage and
        # a WCS object are incomplete.  No way to test the validity of
        # the WCS being copied/created.

        maskedImageMiWcs = self.exposureMiWcs.getMaskedImage()
        miWcsWidth = maskedImageMiWcs.getWidth()
        miWcsHeight = maskedImageMiWcs.getHeight()
        self.assertAlmostEqual(miWcsWidth, self.width)
        self.assertAlmostEqual(miWcsHeight, self.height)

        maskedImageCrWcs = self.exposureCrWcs.getMaskedImage()
        crWcsWidth = maskedImageCrWcs.getWidth()
        crWcsHeight = maskedImageCrWcs.getHeight()
        if crWcsWidth != crWcsHeight != 0:
            self.fail("%s != %s != 0" % (crWcsWidth, crWcsHeight))

        maskedImageCrOnly = self.exposureCrOnly.getMaskedImage()
        crOnlyWidth = maskedImageCrOnly.getWidth()
        crOnlyHeight = maskedImageCrOnly.getHeight()
        if crOnlyWidth != crOnlyHeight != 0:
            self.fail("%s != %s != 0" % (crOnlyWidth, crOnlyHeight))

        # Check Exposure.getWidth() returns the MaskedImage's width
        self.assertEqual(crOnlyWidth, self.exposureCrOnly.getWidth())
        self.assertEqual(crOnlyHeight, self.exposureCrOnly.getHeight())

    def testProperties(self):
        self.assertMaskedImagesEqual(self.exposureMiOnly.maskedImage,
                                     self.exposureMiOnly.getMaskedImage())
        mi2 = afwImage.MaskedImageF(self.exposureMiOnly.getDimensions())
        mi2.image.array[:] = 5.0
        mi2.variance.array[:] = 3.0
        mi2.mask.array[:] = 0x1
        self.exposureMiOnly.maskedImage = mi2
        self.assertMaskedImagesEqual(self.exposureMiOnly.maskedImage, mi2)
        self.assertImagesEqual(self.exposureMiOnly.image,
                               self.exposureMiOnly.maskedImage.image)

        image3 = afwImage.ImageF(self.exposureMiOnly.getDimensions())
        image3.array[:] = 3.0
        self.exposureMiOnly.image = image3
        self.assertImagesEqual(self.exposureMiOnly.image, image3)

        mask3 = afwImage.MaskX(self.exposureMiOnly.getDimensions())
        mask3.array[:] = 0x2
        self.exposureMiOnly.mask = mask3
        self.assertMasksEqual(self.exposureMiOnly.mask, mask3)

        var3 = afwImage.ImageF(self.exposureMiOnly.getDimensions())
        var3.array[:] = 2.0
        self.exposureMiOnly.variance = var3
        self.assertImagesEqual(self.exposureMiOnly.variance, var3)

    def testGetWcs(self):
        """
        Test if a WCS can be obtained from each Exposure created with
        a WCS.

        Test that appropriate exceptions are thrown if a WCS is
        requested from an Exposure that was not created with a WCS.
        Python turns the pex::exceptions in the Exposure and
        MaskedImage classes into IndexErrors.

        The exposureBlank, exposureMiOnly, and exposureCrOnly
        Exposures should throw a lsst::pex::exceptions::NotFound.
        """

        self.assertFalse(self.exposureBlank.getWcs())
        self.assertFalse(self.exposureMiOnly.getWcs())

        # These two should pass
        self.exposureMiWcs.getWcs()
        self.exposureCrWcs.getWcs()

        self.assertFalse(self.exposureCrOnly.getWcs())

    def testExposureInfoConstructor(self):
        """Test the Exposure(maskedImage, exposureInfo) constructor"""
        exposureInfo = afwImage.ExposureInfo()
        exposureInfo.setWcs(self.wcs)
        exposureInfo.setDetector(self.detector)
        gFilter = afwImage.Filter("g")
        exposureInfo.setFilter(gFilter)
        maskedImage = afwImage.MaskedImageF(inFilePathSmall)
        exposure = afwImage.ExposureF(maskedImage, exposureInfo)

        self.assertTrue(exposure.hasWcs())
        self.assertEqual(exposure.getWcs().getPixelOrigin(),
                         self.wcs.getPixelOrigin())
        self.assertEquals(exposure.getDetector().getName(),
                          self.detector.getName())
        self.assertEquals(exposure.getDetector().getSerial(),
                          self.detector.getSerial())
        self.assertEquals(exposure.getFilter(), gFilter)

        self.assertTrue(exposure.getInfo().hasWcs())
        self.assertEqual(exposure.getInfo().getWcs().getPixelOrigin(),
                         self.wcs.getPixelOrigin())
        self.assertEquals(exposure.getInfo().getDetector().getName(),
                          self.detector.getName())
        self.assertEquals(exposure.getInfo().getDetector().getSerial(),
                          self.detector.getSerial())
        self.assertEquals(exposure.getInfo().getFilter(), gFilter)

    def testNullWcs(self):
        """Test that an Exposure constructed with second argument None is usable

        When the exposureInfo constructor was first added, trying to get a WCS
        or other info caused a segfault because the ExposureInfo did not exist.
        """
        maskedImage = self.exposureMiOnly.getMaskedImage()
        exposure = afwImage.ExposureF(maskedImage, None)
        self.assertFalse(exposure.hasWcs())
        self.assertFalse(exposure.hasPsf())

    def testExposureInfoSetNone(self):
        exposureInfo = afwImage.ExposureInfo()
        exposureInfo.setDetector(None)
        exposureInfo.setValidPolygon(None)
        exposureInfo.setPsf(None)
        exposureInfo.setWcs(None)
        exposureInfo.setCalib(None)
        exposureInfo.setCoaddInputs(None)
        exposureInfo.setVisitInfo(None)
        exposureInfo.setApCorrMap(None)

    def testSetExposureInfo(self):
        exposureInfo = afwImage.ExposureInfo()
        exposureInfo.setWcs(self.wcs)
        exposureInfo.setDetector(self.detector)
        gFilter = afwImage.Filter("g")
        exposureInfo.setFilter(gFilter)
        maskedImage = afwImage.MaskedImageF(inFilePathSmall)
        exposure = afwImage.ExposureF(maskedImage)
        self.assertFalse(exposure.hasWcs())

        exposure.setInfo(exposureInfo)

        self.assertTrue(exposure.hasWcs())
        self.assertEqual(exposure.getWcs().getPixelOrigin(),
                         self.wcs.getPixelOrigin())
        self.assertEquals(exposure.getDetector().getName(),
                          self.detector.getName())
        self.assertEquals(exposure.getDetector().getSerial(),
                          self.detector.getSerial())
        self.assertEquals(exposure.getFilter(), gFilter)

    def testVisitInfoFitsPersistence(self):
        """Test saving an exposure to FITS and reading it back in preserves (some) VisitInfo fields"""
        exposureId = 5
        exposureTime = 12.3
        boresightRotAngle = 45.6 * afwGeom.degrees
        weather = afwCoord.Weather(1.1, 2.2, 0.3)
        visitInfo = afwImage.VisitInfo(
            exposureId = exposureId,
            exposureTime = exposureTime,
            boresightRotAngle = boresightRotAngle,
            weather = weather,
        )
        # Calib used to have exposure time and exposure date, so check for lack
        # of interference
        calib = afwImage.Calib(3.4)
        exposureInfo = afwImage.ExposureInfo()
        exposureInfo.setVisitInfo(visitInfo)
        exposureInfo.setCalib(calib)
        exposureInfo.setDetector(self.detector)
        gFilter = afwImage.Filter("g")
        exposureInfo.setFilter(gFilter)
        maskedImage = afwImage.MaskedImageF(inFilePathSmall)
        exposure = afwImage.ExposureF(maskedImage, exposureInfo)
        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            exposure.writeFits(tmpFile)
            rtExposure = afwImage.ExposureF(tmpFile)
        rtVisitInfo = rtExposure.getInfo().getVisitInfo()
        self.assertEqual(rtVisitInfo.getWeather(), weather)
        self.assertEqual(rtExposure.getCalib(), calib)
        self.assertEqual(rtExposure.getFilter(), gFilter)

    def testSetMembers(self):
        """
        Test that the MaskedImage and the WCS of an Exposure can be set.
        """
        exposure = afwImage.ExposureF()

        maskedImage = afwImage.MaskedImageF(inFilePathSmall)
        exposure.setMaskedImage(maskedImage)
        exposure.setWcs(self.wcs)
        exposure.setDetector(self.detector)
        exposure.setFilter(afwImage.Filter("g"))

        self.assertEquals(exposure.getDetector().getName(),
                          self.detector.getName())
        self.assertEquals(exposure.getDetector().getSerial(),
                          self.detector.getSerial())
        self.assertEquals(exposure.getFilter().getName(), "g")

        try:
            exposure.getWcs()
        except pexExcept.Exception as e:
            print("caught expected exception (getWcs): %s" % e)
            pass
        #
        # Test the Calib member.  The Calib tests are in color.py, here we just check that it's in Exposure
        #
        calib = exposure.getCalib()
        fluxMag0, FluxMag0Err = 1.1e12, 2.2e10
        calib.setFluxMag0(fluxMag0, FluxMag0Err)
        self.assertEqual(exposure.getCalib().getFluxMag0(),
                         (fluxMag0, FluxMag0Err))
        #
        # now check that we can set Calib
        #
        calib = afwImage.Calib()
        fluxMag0_2, fluxMag0Err_2 = (511.1, 44.4)
        calib.setFluxMag0(fluxMag0_2, fluxMag0Err_2)

        exposure.setCalib(calib)

        self.assertEqual(exposure.getCalib().getFluxMag0(),
                         (fluxMag0_2, fluxMag0Err_2))
        #
        # Psfs next
        #
        self.assertFalse(exposure.hasPsf())
        exposure.setPsf(self.psf)
        self.assertTrue(exposure.hasPsf())

        exposure.setPsf(DummyPsf(1.0))  # we can reset the Psf

        # Test that we can set the MaskedImage and WCS of an Exposure
        # that already has both
        self.exposureMiWcs.setMaskedImage(maskedImage)
        exposure.setWcs(self.wcs)

    def testHasWcs(self):
        """
        Test if an Exposure has a WCS or not.
        """
        self.assertFalse(self.exposureBlank.hasWcs())

        self.assertFalse(self.exposureMiOnly.hasWcs())
        self.assertTrue(self.exposureMiWcs.hasWcs())
        self.assertTrue(self.exposureCrWcs.hasWcs())
        self.assertFalse(self.exposureCrOnly.hasWcs())

    def testGetSubExposure(self):
        """
        Test that a subExposure of the original Exposure can be obtained.

        The MaskedImage class should throw a
        lsst::pex::exceptions::InvalidParameter if the requested
        subRegion is not fully contained within the original
        MaskedImage.

        """
        #
        # This subExposure is valid
        #
        subBBox = afwGeom.Box2I(afwGeom.Point2I(40, 50),
                                afwGeom.Extent2I(10, 10))
        subExposure = self.exposureCrWcs.Factory(
            self.exposureCrWcs, subBBox, afwImage.LOCAL)

        self.checkWcs(self.exposureCrWcs, subExposure)

        # this subRegion is not valid and should trigger an exception
        # from the MaskedImage class and should trigger an exception
        # from the WCS class for the MaskedImage 871034p_1_MI.

        subRegion3 = afwGeom.Box2I(afwGeom.Point2I(100, 100),
                                   afwGeom.Extent2I(10, 10))

        def getSubRegion():
            self.exposureCrWcs.Factory(
                self.exposureCrWcs, subRegion3, afwImage.LOCAL)

        self.assertRaises(pexExcept.LengthError, getSubRegion)

        # this subRegion is not valid and should trigger an exception
        # from the MaskedImage class only for the MaskedImage small_MI.
        # small_MI (cols, rows) = (256, 256)

        subRegion4 = afwGeom.Box2I(afwGeom.Point2I(250, 250),
                                   afwGeom.Extent2I(10, 10))

        def getSubRegion():
            self.exposureCrWcs.Factory(
                self.exposureCrWcs, subRegion4, afwImage.LOCAL)

        self.assertRaises(pexExcept.LengthError, getSubRegion)

        # check the sub- and parent- exposures are using the same Wcs
        # transformation
        subBBox = afwGeom.Box2I(afwGeom.Point2I(40, 50),
                                afwGeom.Extent2I(10, 10))
        subExposure = self.exposureCrWcs.Factory(
            self.exposureCrWcs, subBBox, afwImage.LOCAL)
        parentPos = self.exposureCrWcs.getWcs().pixelToSky(0, 0)

        parentPos = parentPos.getPosition()

        subExpPos = subExposure.getWcs().pixelToSky(0, 0).getPosition()

        for i in range(2):
            self.assertAlmostEqual(
                parentPos[i], subExpPos[i], 9, "Wcs in sub image has changed")

    def testReadWriteFits(self):
        """Test readFits and writeFits.
        """
        # This should pass without an exception
        mainExposure = afwImage.ExposureF(inFilePathSmall)
        mainExposure.setDetector(self.detector)

        subBBox = afwGeom.Box2I(afwGeom.Point2I(10, 10),
                                afwGeom.Extent2I(40, 50))
        subExposure = mainExposure.Factory(
            mainExposure, subBBox, afwImage.LOCAL)
        self.checkWcs(mainExposure, subExposure)
        det = subExposure.getDetector()
        self.assertTrue(det)

        subExposure = afwImage.ExposureF(
            inFilePathSmall, subBBox, afwImage.LOCAL)

        self.checkWcs(mainExposure, subExposure)

        # This should throw an exception
        def getExposure():
            afwImage.ExposureF(inFilePathSmallImage)

        self.assertRaises(lsst.afw.fits.FitsError, getExposure)

        mainExposure.setPsf(self.psf)

        # Make sure we can write without an exception
        fluxMag0, fluxMag0Err = 1e12, 1e10
        mainExposure.getCalib().setFluxMag0(fluxMag0, fluxMag0Err)

        # Check scaling of Calib
        scale = 2.0
        calib = mainExposure.getCalib()
        self.assertEqual((fluxMag0, fluxMag0Err), calib.getFluxMag0())
        self.assertEqual((fluxMag0, fluxMag0Err),
                         mainExposure.getCalib().getFluxMag0())
        calib *= scale
        self.assertEqual((fluxMag0*scale, fluxMag0Err*scale),
                         calib.getFluxMag0())
        self.assertEqual((fluxMag0*scale, fluxMag0Err*scale),
                         mainExposure.getCalib().getFluxMag0())
        calib /= scale
        self.assertEqual((fluxMag0, fluxMag0Err), calib.getFluxMag0())
        self.assertEqual((fluxMag0, fluxMag0Err),
                         mainExposure.getCalib().getFluxMag0())

        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            mainExposure.writeFits(tmpFile)

            readExposure = type(mainExposure)(tmpFile)

            #
            # Check the round-tripping
            #
            self.assertEqual(mainExposure.getFilter().getName(),
                             readExposure.getFilter().getName())

            self.assertEqual((fluxMag0, fluxMag0Err),
                             readExposure.getCalib().getFluxMag0())

            psf = readExposure.getPsf()
            self.assertIsNotNone(psf)
            self.assertEqual(psf.getValue(), self.psf.getValue())

    def checkWcs(self, parentExposure, subExposure):
        """Compare WCS at corner points of a sub-exposure and its parent exposure
           By using the function indexToPosition, we should be able to convert the indices
           (of the four corners (of the sub-exposure)) to positions and use the wcs
           to get the same sky coordinates for each.
        """
        subMI = subExposure.getMaskedImage()
        subDim = subMI.getDimensions()

        # Note: pixel positions must be computed relative to XY0 when working
        # with WCS
        mainWcs = parentExposure.getWcs()
        subWcs = subExposure.getWcs()

        for xSubInd in (0, subDim.getX()-1):
            for ySubInd in (0, subDim.getY()-1):
                mainWcs.pixelToSky(
                    afwImage.indexToPosition(xSubInd),
                    afwImage.indexToPosition(ySubInd),
                )
                subWcs.pixelToSky(
                    afwImage.indexToPosition(xSubInd),
                    afwImage.indexToPosition(ySubInd),
                )

    def cmpExposure(self, e1, e2):
        self.assertEqual(e1.getDetector().getName(),
                         e2.getDetector().getName())
        self.assertEqual(e1.getDetector().getSerial(),
                         e2.getDetector().getSerial())
        self.assertEqual(e1.getFilter().getName(), e2.getFilter().getName())
        xy = afwGeom.Point2D(0, 0)
        self.assertEqual(e1.getWcs().pixelToSky(xy)[0],
                         e2.getWcs().pixelToSky(xy)[0])
        self.assertEqual(e1.getCalib(), e2.getCalib())
        # check PSF identity
        if not e1.getPsf():
            self.assertFalse(e2.getPsf())
        else:
            self.assertEqual(e1.getPsf().getValue(), e2.getPsf().getValue())

    def testCopyExposure(self):
        """Copy an Exposure (maybe changing type)"""

        exposureU = afwImage.ExposureU(inFilePathSmall)
        exposureU.setWcs(self.wcs)
        exposureU.setDetector(self.detector)
        exposureU.setFilter(afwImage.Filter("g"))
        exposureU.setPsf(DummyPsf(4.0))

        exposureF = exposureU.convertF()
        self.cmpExposure(exposureF, exposureU)

        nexp = exposureF.Factory(exposureF, False)
        self.cmpExposure(exposureF, nexp)

        # Ensure that the copy was deep.
        # (actually this test is invalid since getDetector() returns a shared_ptr)
        # cen0 = exposureU.getDetector().getCenterPixel()
        # x0,y0 = cen0
        # det = exposureF.getDetector()
        # det.setCenterPixel(afwGeom.Point2D(999.0, 437.8))
        # self.assertEqual(exposureU.getDetector().getCenterPixel()[0], x0)
        # self.assertEqual(exposureU.getDetector().getCenterPixel()[1], y0)

    def testDeepCopyData(self):
        """Make sure a deep copy of an Exposure has its own data (ticket #2625)
        """
        exp = afwImage.ExposureF(6, 7)
        mi = exp.getMaskedImage()
        mi.getImage().set(100)
        mi.getMask().set(5)
        mi.getVariance().set(200)

        expCopy = exp.clone()
        miCopy = expCopy.getMaskedImage()
        miCopy.getImage().set(-50)
        miCopy.getMask().set(2)
        miCopy.getVariance().set(175)

        self.assertFloatsAlmostEqual(miCopy.getImage().getArray(), -50)
        self.assertTrue(np.all(miCopy.getMask().getArray() == 2))
        self.assertFloatsAlmostEqual(miCopy.getVariance().getArray(), 175)

        self.assertFloatsAlmostEqual(mi.getImage().getArray(), 100)
        self.assertTrue(np.all(mi.getMask().getArray() == 5))
        self.assertFloatsAlmostEqual(mi.getVariance().getArray(), 200)

    def testDeepCopySubData(self):
        """Make sure a deep copy of a subregion of an Exposure has its own data (ticket #2625)
        """
        exp = afwImage.ExposureF(6, 7)
        mi = exp.getMaskedImage()
        mi.getImage().set(100)
        mi.getMask().set(5)
        mi.getVariance().set(200)

        bbox = afwGeom.Box2I(afwGeom.Point2I(1, 0), afwGeom.Extent2I(5, 4))
        expCopy = exp.Factory(exp, bbox, afwImage.PARENT, True)
        miCopy = expCopy.getMaskedImage()
        miCopy.getImage().set(-50)
        miCopy.getMask().set(2)
        miCopy.getVariance().set(175)

        self.assertFloatsAlmostEqual(miCopy.getImage().getArray(), -50)
        self.assertTrue(np.all(miCopy.getMask().getArray() == 2))
        self.assertFloatsAlmostEqual(miCopy.getVariance().getArray(), 175)

        self.assertFloatsAlmostEqual(mi.getImage().getArray(), 100)
        self.assertTrue(np.all(mi.getMask().getArray() == 5))
        self.assertFloatsAlmostEqual(mi.getVariance().getArray(), 200)

    def testDeepCopyMetadata(self):
        """Make sure a deep copy of an Exposure has a deep copy of metadata (ticket #2568)
        """
        exp = afwImage.ExposureF(10, 10)
        expMeta = exp.getMetadata()
        expMeta.set("foo", 5)
        expCopy = exp.clone()
        expCopyMeta = expCopy.getMetadata()
        expCopyMeta.set("foo", 6)
        self.assertEqual(expCopyMeta.get("foo"), 6)
        # this will fail if the bug is present
        self.assertEqual(expMeta.get("foo"), 5)

    def testDeepCopySubMetadata(self):
        """Make sure a deep copy of a subregion of an Exposure has a deep copy of metadata (ticket #2568)
        """
        exp = afwImage.ExposureF(10, 10)
        expMeta = exp.getMetadata()
        expMeta.set("foo", 5)
        bbox = afwGeom.Box2I(afwGeom.Point2I(1, 0), afwGeom.Extent2I(5, 5))
        expCopy = exp.Factory(exp, bbox, afwImage.PARENT, True)
        expCopyMeta = expCopy.getMetadata()
        expCopyMeta.set("foo", 6)
        self.assertEqual(expCopyMeta.get("foo"), 6)
        # this will fail if the bug is present
        self.assertEqual(expMeta.get("foo"), 5)

    def testMakeExposureLeaks(self):
        """Test for memory leaks in makeExposure (the test is in lsst.utils.tests.MemoryTestCase)"""
        afwImage.makeMaskedImage(afwImage.ImageU(afwGeom.Extent2I(10, 20)))
        afwImage.makeExposure(afwImage.makeMaskedImage(
            afwImage.ImageU(afwGeom.Extent2I(10, 20))))

    def testImageSlices(self):
        """Test image slicing, which generate sub-images using Box2I under the covers"""
        exp = afwImage.ExposureF(10, 20)
        mi = exp.getMaskedImage()
        mi[9, 19] = 10
        # N.b. Exposures don't support setting/getting the pixels so can't
        # replicate e.g. Image's slice tests
        sexp = exp[1:4, 6:10]
        self.assertEqual(sexp.getDimensions(), afwGeom.ExtentI(3, 4))
        sexp = exp[..., -3:]
        self.assertEqual(sexp.getDimensions(),
                         afwGeom.ExtentI(exp.getWidth(), 3))
        self.assertEqual(sexp.getMaskedImage().get(sexp.getWidth() - 1, sexp.getHeight() - 1),
                         exp.getMaskedImage().get(exp.getWidth() - 1, exp.getHeight() - 1))

    def testConversionToScalar(self):
        """Test that even 1-pixel Exposures can't be converted to scalars"""
        im = afwImage.ExposureF(10, 20)

        # only single pixel images may be converted
        self.assertRaises(TypeError, float, im)
        # actually, can't convert (img, msk, var) to scalar
        self.assertRaises(TypeError, float, im[0, 0])

    def testReadMetadata(self):
        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            self.exposureCrWcs.getMetadata().set("FRAZZLE", True)
            # This will write the main metadata (inc. FRAZZLE) to the primary HDU, and the
            # WCS to subsequent HDUs, along with INHERIT=T.
            self.exposureCrWcs.writeFits(tmpFile)
            # This should read the first non-empty HDU (i.e. it skips the primary), but
            # goes back and reads it if it finds INHERIT=T.  That should let us read
            # frazzle and the Wcs from the PropertySet returned by
            # readMetadata.
            md = afwImage.readMetadata(tmpFile)
            wcs = afwImage.makeWcs(md, True)
            self.assertEqual(wcs.getPixelOrigin(), self.wcs.getPixelOrigin())
            self.assertEqual(wcs.getSkyOrigin(), self.wcs.getSkyOrigin())
            self.assertTrue(
                np.all(wcs.getCDMatrix() == self.wcs.getCDMatrix()))
            frazzle = md.get("FRAZZLE")
            self.assertTrue(frazzle)

    def testArchiveKeys(self):
        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            exposure1 = afwImage.ExposureF(100, 100, self.wcs)
            exposure1.setPsf(self.psf)
            exposure1.writeFits(tmpFile)
            exposure2 = afwImage.ExposureF(tmpFile)
            self.assertFalse(exposure2.getMetadata().exists("AR_ID"))
            self.assertFalse(exposure2.getMetadata().exists("PSF_ID"))
            self.assertFalse(exposure2.getMetadata().exists("WCS_ID"))

    def testTicket2861(self):
        with lsst.utils.tests.getTempFilePath(".fits") as tmpFile:
            exposure1 = afwImage.ExposureF(100, 100, self.wcs)
            exposure1.setPsf(self.psf)
            schema = afwTable.ExposureTable.makeMinimalSchema()
            coaddInputs = afwImage.CoaddInputs(schema, schema)
            exposure1.getInfo().setCoaddInputs(coaddInputs)
            exposure2 = afwImage.ExposureF(exposure1, True)
            self.assertIsNotNone(exposure2.getInfo().getCoaddInputs())
            exposure2.writeFits(tmpFile)
            exposure3 = afwImage.ExposureF(tmpFile)
            self.assertIsNotNone(exposure3.getInfo().getCoaddInputs())


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
