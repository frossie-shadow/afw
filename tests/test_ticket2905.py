#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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
import os.path
import unittest

import lsst.afw.image as afwImage
import lsst.utils.tests

testPath = os.path.abspath(os.path.dirname(__file__))


class Ticket2905Test(unittest.TestCase):
    """Test reading a FITS header that contains:

    INR-STR =                2E-05
    """

    def test(self):
        path = os.path.join(testPath, "data", "ticket2905.fits")
        md = afwImage.readMetadata(path)
        value = md.get("INR-STR")
        self.assertEqual(type(value), float)
        self.assertEqual(value, 2.0e-5)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
