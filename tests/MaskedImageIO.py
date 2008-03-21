"""
Tests for MaskedImages

Run with:
   python MaskedImageIO.py
or
   python
   >>> import MaskedImageIO; MaskedImageIO.run()
"""

import pdb                              # we may want to say pdb.set_trace()
import os
import eups
import unittest
import lsst.mwi.tests as tests
import lsst.fw.Core.fwLib as fw
import lsst.fw.Display.ds9 as ds9
import fwTests
from lsst.mwi.exceptions import *

dataDir = eups.productDir("fwData")
if not dataDir:
    raise RuntimeError("You must set up fwData to run these tests")

try:
    type(display)
except NameError:
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class MaskedImageTestCase(unittest.TestCase):
    """A test case for MaskedImage"""
    def setUp(self):
        self.mi = fw.MaskedImageF()

        mask = self.mi.getMask()

        # Set a (non-standard) initial Mask plane definition
        #
        # Ideally we'd use the standard dictionary and a non-standard file, but
        # a standard file's what we have
        #
        mask.clearMaskPlaneDict()
        for p in ("ZERO", "BAD", "SAT", "INTRP", "CR", "EDGE"):
            mask.addMaskPlane(p)

        self.file = os.path.join(dataDir, "CFHT", "D4", "cal-53535-i-797722_1")

    def tearDown(self):
        del self.mi

    def testFitsRead(self):
        """Check if we read MaskedImages"""

        image = self.mi.getImage()
        mask = self.mi.getMask()

        self.mi.readFits(self.file)

        #mask.printMaskPlanes()

        if display:
            ds9.mtv(self.mi)

        self.assertEqual(image.getVal(32, 1), 3728);
        self.assertEqual(mask.getVal(0,0), 2); # BAD
            
    def testFitsReadNoConform(self):
        """Check if we read MaskedImages and make them conform to Mask's plane dictionary"""

        image = self.mi.getImage()
        mask = self.mi.getMask()

        self.mi.readFits(self.file, False)

        self.assertEqual(image.getVal(32, 1), 3728);
        self.assertEqual(mask.getVal(0,0), 2); # i.e. shifted 1 place to the right

        self.assertEqual(mask.getMaskPlane("CR"), 4, "Plane CR has value defined by Mask")

    def testFitsReadConform(self):
        """Check if we read MaskedImages and make them replace Mask's plane dictionary"""

        image = self.mi.getImage()
        mask = self.mi.getMask()

        self.mi.readFits(self.file, True)

        self.assertEqual(image.getVal(32, 1), 3728);
        self.assertEqual(mask.getVal(0,0), 1); # i.e. not shifted 1 place to the right

        self.assertEqual(mask.getMaskPlane("CR"), 3, "Plane CR has value specified in FITS file")

    def testFitsReadNoConform2(self):
        """Check that reading a mask doesn't invalidate the plane dictionary"""

        mask = self.mi.getMask()

        self.mi.readFits(self.file)     # just to get the size
        testMask = fw.MaskU(mask.getCols(), mask.getRows())
        self.mi.readFits(self.file, False)

        mask |= testMask

    def testFitsReadConform2(self):
        """Check that conforming a mask invalidates the plane dictionary"""

        mask = self.mi.getMask()

        self.mi.readFits(self.file)     # just to get the size
        testMask = fw.MaskU(mask.getCols(), mask.getRows())
        self.mi.readFits(self.file, True)

        def tst(mask=mask):
            mask |= testMask

        self.assertRaises(LsstRuntime, tst)
        mask.this.acquire() # Work around bug in swig "mask |= mask;" leaks when |= throws

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""

    tests.init()

    suites = []
    suites += unittest.makeSuite(MaskedImageTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    tests.run(suite(), exit)        # mwi 2.0

if __name__ == "__main__":
    run(True)
