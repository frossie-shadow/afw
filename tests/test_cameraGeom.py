
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

from __future__ import absolute_import, division, print_function
import unittest
import os

from builtins import next
from builtins import zip
from builtins import range
import numpy as np

import lsst.utils.tests
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9
from lsst.afw.cameraGeom import PIXELS, FIELD_ANGLE, FOCAL_PLANE, CameraSys, CameraSysPrefix, \
    CameraPoint, Camera, Detector, assembleAmplifierImage, assembleAmplifierRawImage
import lsst.afw.cameraGeom.testUtils as testUtils
import lsst.afw.cameraGeom.utils as cameraGeomUtils

try:
    type(display)
except NameError:
    display = False


testPath = os.path.abspath(os.path.dirname(__file__))


class CameraGeomTestCase(unittest.TestCase):
    """A test case for camera geometry"""

    def setUp(self):
        self.lsstCamWrapper = testUtils.CameraWrapper(isLsstLike=True)
        self.scCamWrapper = testUtils.CameraWrapper(isLsstLike=False)
        self.cameraList = (self.lsstCamWrapper, self.scCamWrapper)
        self.assemblyList = {}
        self.assemblyList[self.lsstCamWrapper.camera.getName()] =\
            [afwImage.ImageU(os.path.join(testPath, 'test_amp.fits.gz'))
             for i in range(8)]
        self.assemblyList[self.scCamWrapper.camera.getName()] =\
            [afwImage.ImageU(os.path.join(testPath, 'test.fits.gz'))]

    def tearDown(self):
        del self.lsstCamWrapper
        del self.scCamWrapper
        del self.cameraList
        del self.assemblyList

    def testConstructor(self):
        for cw in self.cameraList:
            self.assertIsInstance(cw.camera, Camera)
            self.assertEqual(cw.nDetectors, len(cw.camera))
            self.assertEqual(cw.nDetectors, len(cw.ampInfoDict))
            self.assertEqual(sorted(cw.detectorNameList),
                             sorted(cw.camera.getNameIter()))
            self.assertEqual(sorted(cw.detectorIdList),
                             sorted(cw.camera.getIdIter()))
            for det in cw.camera:
                self.assertIsInstance(det, Detector)
                self.assertEqual(
                    cw.ampInfoDict[det.getName()]['namps'], len(det))

    def testMakeCameraPoint(self):
        point = afwGeom.Point2D(0, 0)
        for cw in self.cameraList:
            camera = cw.camera
            for coordSys in (FIELD_ANGLE, FOCAL_PLANE):
                pt1 = afwGeom.Point2D(0.1, 0.3)
                pt2 = afwGeom.Point2D(0., 0.)
                pt3 = afwGeom.Point2D(-0.2, 0.2)
                pt4 = afwGeom.Point2D(0.02, -0.2)
                pt5 = afwGeom.Point2D(-0.2, -0.03)
                for pt in (pt1, pt2, pt3, pt4, pt5):
                    cp = camera.makeCameraPoint(pt, coordSys)
                    self.assertEquals(cp.getPoint(), pt)
                    self.assertEquals(
                        cp.getCameraSys().getSysName(), coordSys.getSysName())

                    # test == and !=
                    cp2 = camera.makeCameraPoint(pt, coordSys)
                    self.assertTrue(cp == cp2)
                    self.assertFalse(cp != cp2)

            det = camera[next(camera.getNameIter())]
            cp = camera.makeCameraPoint(point, FOCAL_PLANE)
            self.checkCamPoint(cp, point, FOCAL_PLANE)
            cp = camera.makeCameraPoint(point, det.makeCameraSys(PIXELS))
            self.checkCamPoint(cp, point, det.makeCameraSys(PIXELS))
            # non-existant camera sys in makeCameraPoint
            self.assertRaises(
                RuntimeError, camera.makeCameraPoint, point, CameraSys('abcd'))
            # CameraSysPrefix camera sys in makeCameraPoint
            self.assertRaises(TypeError, camera.makeCameraPoint, point, PIXELS)

    def testCameraSysRepr(self):
        """Test CameraSys.__repr__ and CameraSysPrefix.__repr__
        """
        for sysName in ("FocalPlane", "FieldAngle", "Pixels", "foo"):
            cameraSys = CameraSys(sysName)
            predRepr = "CameraSys(%s)" % (sysName)
            self.assertEqual(repr(cameraSys), predRepr)

            cameraSysPrefix = CameraSysPrefix(sysName)
            predCSPRepr = "CameraSysPrefix(%s)" % (sysName)
            self.assertEqual(repr(cameraSysPrefix), predCSPRepr)
            for detectorName in ("Detector 1", "bar"):
                cameraSys2 = CameraSys(sysName, detectorName)
                predRepr2 = "CameraSys(%s, %s)" % (sysName, detectorName)
                self.assertEqual(repr(cameraSys2), predRepr2)

    def testCameraPointRepr(self):
        """Test CameraPoint.__repr__
        """
        point = afwGeom.Point2D(1.5, -23.4)
        cameraSys = FOCAL_PLANE
        cameraPoint = CameraPoint(point, cameraSys)
        predRepr = "CameraPoint(%s, %s)" % (point, cameraSys)
        self.assertEqual(repr(cameraPoint), predRepr)

    def testAccessor(self):
        for cw in self.cameraList:
            camera = cw.camera
            for name in cw.detectorNameList:
                self.assertIsInstance(camera[name], Detector)
            for detId in cw.detectorIdList:
                self.assertIsInstance(camera[detId], Detector)

    def testTransformSlalib(self):
        """Test Camera.transform against data computed using SLALIB

        These test data come from SLALIB using SLA_PCD with 0.925 and
        a plate scale of 20 arcsec/mm
        """
        testData = [(-1.84000000, 1.04000000, -331.61689069, 187.43563387),
                    (-1.64000000, 0.12000000, -295.42491556, 21.61645724),
                    (-1.44000000, -0.80000000, -259.39818797, -144.11010443),
                    (-1.24000000, -1.72000000, -223.48275934, -309.99221457),
                    (-1.08000000, 1.36000000, -194.56520533, 245.00803635),
                    (-0.88000000, 0.44000000, -158.44320430, 79.22160215),
                    (-0.68000000, -0.48000000, -122.42389383, -86.41686623),
                    (-0.48000000, -1.40000000, -86.45332534, -252.15553224),
                    (-0.32000000, 1.68000000, -57.64746955, 302.64921514),
                    (-0.12000000, 0.76000000, -21.60360306, 136.82281940),
                    (0.08000000, -0.16000000, 14.40012984, -28.80025968),
                    (0.28000000, -1.08000000, 50.41767773, -194.46818554),
                    (0.48000000, -2.00000000, 86.50298919, -360.42912163),
                    (0.64000000, 1.08000000, 115.25115701, 194.48632746),
                    (0.84000000, 0.16000000, 151.23115189, 28.80593369),
                    (1.04000000, -0.76000000, 187.28751874, -136.86395600),
                    (1.24000000, -1.68000000, 223.47420612, -302.77150507),
                    (1.40000000, 1.40000000, 252.27834478, 252.27834478),
                    (1.60000000, 0.48000000, 288.22644118, 86.46793236),
                    (1.80000000, -0.44000000, 324.31346653, -79.27662515), ]

        for cw in self.cameraList:
            camera = cw.camera
            for point in testData:
                fpGivenPos = afwGeom.Point2D(point[2], point[3])
                fpGivenCP = camera.makeCameraPoint(fpGivenPos, FOCAL_PLANE)
                fieldGivenPos = afwGeom.Point2D(
                    afwGeom.degToRad(point[0]), afwGeom.degToRad(point[1]))
                fieldGivenCP = camera.makeCameraPoint(fieldGivenPos, FIELD_ANGLE)

                fpComputedCP = camera.transform(fieldGivenCP, FOCAL_PLANE)
                self.assertCamPointAlmostEquals(fpComputedCP, fpGivenCP)

                fieldComputedCP = camera.transform(fpGivenCP, FIELD_ANGLE)
                self.assertCamPointAlmostEquals(fieldComputedCP, fieldGivenCP)

    def testTransformDet(self):
        """Test Camera.transform with detector-based coordinate systems (PIXELS)
        """
        for cw in self.cameraList:
            numOffUsable = 0  # number of points off one detector but on another
            camera = cw.camera
            detNameList = list(camera.getNameIter())
            for detName in detNameList:
                det = camera[detName]

                # test transforms using a point on the detector
                pixCP = det.makeCameraPoint(afwGeom.Point2D(10, 10), PIXELS)
                fpCP = camera.transform(pixCP, FOCAL_PLANE)
                fieldCP = camera.transform(pixCP, FIELD_ANGLE)

                fieldCP2 = camera.transform(fpCP, FIELD_ANGLE)
                self.assertCamPointAlmostEquals(fieldCP, fieldCP2)

                for intermedCP in (pixCP, fpCP, fieldCP):
                    pixRoundTripCP = camera.transform(
                        intermedCP, det.makeCameraSys(PIXELS))
                    self.assertCamPointAlmostEquals(pixCP, pixRoundTripCP)

                    pixFindRoundTripCP = camera.transform(intermedCP, PIXELS)
                    self.assertCamPointAlmostEquals(pixCP, pixFindRoundTripCP)

                # test transforms using a point off the detector
                pixOffDetCP = det.makeCameraPoint(
                    afwGeom.Point2D(0, -10), PIXELS)
                pixOffDetRoundTripCP = camera.transform(
                    pixOffDetCP, det.makeCameraSys(PIXELS))
                self.assertCamPointAlmostEquals(
                    pixOffDetCP, pixOffDetRoundTripCP)

                # the point off the detector MAY be on another detector
                # (depending if the detector has neighbor on the correct edge)
                detList = camera.findDetectors(pixOffDetCP)
                if len(detList) == 1:
                    numOffUsable += 1
                    pixFindOffCP = camera.transform(pixOffDetCP, PIXELS)
                    self.assertNotEqual(
                        pixCP.getCameraSys(), pixFindOffCP.getCameraSys())

                    # convert point on other detector to pixels on the main detector
                    # the result should not be on the main detector
                    pixToPixCP = camera.transform(
                        pixFindOffCP, det.makeCameraSys(PIXELS))
                    self.assertFalse(afwGeom.Box2D(
                        det.getBBox()).contains(pixToPixCP.getPoint()))
            self.assertGreater(numOffUsable, 0)
            print("numOffUsable=", numOffUsable)

    def testFindDetectors(self):
        for cw in self.cameraList:
            detPointsList = []
            for det in cw.camera:
                # This currently assumes there is only one detector at the center
                # position of any detector.  That is not enforced and multiple detectors
                # at a given FIELD_ANGLE position is supported.  Change this if the default
                # camera changes.
                cp = det.getCenter(FOCAL_PLANE)
                detPointsList.append(cp.getPoint())
                detList = cw.camera.findDetectors(cp)
                self.assertEquals(len(detList), 1)
                self.assertEquals(det.getName(), detList[0].getName())
            detList = cw.camera.findDetectorsList(detPointsList, FOCAL_PLANE)
            self.assertEquals(len(cw.camera), len(detList))
            for dets in detList:
                self.assertEquals(len(dets), 1)

    def testFpBbox(self):
        for cw in self.cameraList:
            camera = cw.camera
            bbox = afwGeom.Box2D()
            for name in cw.detectorNameList:
                for corner in camera[name].getCorners(FOCAL_PLANE):
                    bbox.include(corner)
            self.assertEqual(bbox.getMin(), camera.getFpBBox().getMin())
            self.assertEqual(bbox.getMax(), camera.getFpBBox().getMax())

    def testLinearity(self):
        """Test if we can set/get Linearity parameters"""
        for cw in self.cameraList:
            camera = cw.camera
            for det in camera:
                for amp in det:
                    self.assertEquals(cw.ampInfoDict[det.getName()]['linInfo'][amp.getName()]['linthresh'],
                                      amp.get('linearityThreshold'))
                    self.assertEquals(cw.ampInfoDict[det.getName()]['linInfo'][amp.getName()]['linmax'],
                                      amp.get('linearityMaximum'))
                    self.assertEquals(cw.ampInfoDict[det.getName()]['linInfo'][amp.getName()]['linunits'],
                                      amp.get('linearityUnits'))
                    self.assertEquals(cw.ampInfoDict[det.getName()]['linInfo'][amp.getName()]['lintype'],
                                      amp.getLinearityType())
                    for c1, c2 in zip(cw.ampInfoDict[det.getName()]['linInfo'][amp.getName()]['lincoeffs'],
                                      amp.getLinearityCoeffs()):
                        if np.isfinite(c1) and np.isfinite(c2):
                            self.assertEquals(c1, c2)

    def testAssembly(self):
        ccdNames = ('R:0,0 S:1,0', 'R:0,0 S:0,1')
        compMap = {True: afwImage.ImageU(os.path.join(testPath, 'test_comp_trimmed.fits.gz')),
                   False: afwImage.ImageU(os.path.join(testPath, 'test_comp.fits.gz'))}
        for cw in self.cameraList:
            camera = cw.camera
            imList = self.assemblyList[camera.getName()]
            for ccdName in ccdNames:
                for trim, assemble in ((False, assembleAmplifierRawImage), (True, assembleAmplifierImage)):
                    det = camera[ccdName]
                    if not trim:
                        outBbox = cameraGeomUtils.calcRawCcdBBox(det)
                    else:
                        outBbox = det.getBBox()
                    outImage = afwImage.ImageU(outBbox)
                    if len(imList) == 1:
                        for amp in det:
                            assemble(outImage, imList[0], amp)
                    else:
                        for amp, im in zip(det, imList):
                            assemble(outImage, im, amp)
                    self.assertListEqual(outImage.getArray().flatten().tolist(),
                                         compMap[trim].getArray().flatten().tolist())

    @unittest.skipIf(not display, "display variable not set; skipping cameraGeomUtils test")
    def testCameraGeomUtils(self):
        for cw in self.cameraList:
            camera = cw.camera
            cameraGeomUtils.showCamera(
                camera, referenceDetectorName=camera[0].getName())
            ds9.incrDefaultFrame()
            for det in (camera[0], camera[1]):
                cameraGeomUtils.showCcd(
                    det, isTrimmed=True, inCameraCoords=False)
                ds9.incrDefaultFrame()
                cameraGeomUtils.showCcd(
                    det, isTrimmed=True, inCameraCoords=True)
                ds9.incrDefaultFrame()
                cameraGeomUtils.showCcd(
                    det, isTrimmed=False, inCameraCoords=False)
                ds9.incrDefaultFrame()
                cameraGeomUtils.showCcd(
                    det, isTrimmed=False, inCameraCoords=True)
                ds9.incrDefaultFrame()
                for amp in det:
                    cameraGeomUtils.showAmp(amp)
                    ds9.incrDefaultFrame()

    def testCameraRaises(self):
        for cw in self.cameraList:
            camera = cw.camera
            cp = camera.makeCameraPoint(afwGeom.Point2D(1e6, 1e6), FOCAL_PLANE)
            # Way off the focal plane
            self.assertRaises(RuntimeError, camera.transform, cp, PIXELS)
            # non-existant destination camera system
            cp = camera.makeCameraPoint(afwGeom.Point2D(0, 0), FOCAL_PLANE)
            self.assertRaises(RuntimeError, camera.transform,
                              cp, CameraSys('abcd'))

    def checkCamPoint(self, cp, testPt, testSys):
        """Assert that a CameraPoint contains the specified Point2D and CameraSys"""
        self.assertEquals(cp.getCameraSys(), testSys)
        self.assertEquals(cp.getPoint(), testPt)

    def assertCamPointAlmostEquals(self, cp1, cp2, ndig=6):
        """Assert that two CameraPoints are nearly equal
        """
        self.assertEquals(cp1.getCameraSys(), cp2.getCameraSys())
        for i in range(2):
            self.assertAlmostEquals(cp1.getPoint()[i], cp2.getPoint()[i], 6)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
