# This file is part of obs_lsst.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""Test support classes for obs_lsst"""
__all__ = ("ObsLsstButlerTests", "ObsLsstObsBaseOverrides")

import os.path
import unittest
import lsst.utils.tests
import lsst.obs.base.tests
from lsst.utils import getPackageDir

# Define the data location relative to this package
DATAROOT = os.path.join(getPackageDir("obs_lsst"), "data", "input")


class ObsLsstButlerTests(lsst.utils.tests.TestCase):
    """Base class shared by all tests of the butler and mapper.

    This class can not inherit from `~lsst.obs.base.tests.ObsTests` since
    that will trigger tests in this class directly that will fail.

    This class defines a butler and a mapper for each test subclass.
    They are stored in the ``_mapper`` and ``_butler`` class attributes
    to distinguish them from the ``mapper`` and ``butler`` instance
    attributes used by `~lsst.obs.base.tests.ObsTests`.
    """

    instrumentDir = "TBD"  # Override in subclass
    """Name of instrument directory within data/input."""

    _mapper = None
    _butler = None

    @classmethod
    def tearDownClass(cls):
        del cls._mapper
        del cls._butler

    @classmethod
    def setUpClass(cls):
        cls.data_dir = os.path.normpath(os.path.join(DATAROOT, cls.instrumentDir))
        # Protection against the base class values being used
        if not os.path.exists(cls.data_dir):
            raise unittest.SkipTest(f"Data directory {cls.data_dir} does not exist.")

        cls._butler = lsst.daf.persistence.Butler(root=cls.data_dir)
        mapper_class = cls._butler.getMapperClass(root=cls.data_dir)
        cls._mapper = mapper_class(root=cls.data_dir)


class ObsLsstObsBaseOverrides(lsst.obs.base.tests.ObsTests):
    """Specialist butler tests for obs_lsst."""

    @unittest.skip("raw_sub not supported bt obs_lsst")
    def test_raw_sub_bbox(self):
        return

    def testRawVisitInfo(self):
        visitInfo = self.butler.get("raw_visitInfo", self.dataIds["raw"])
        self.assertIsInstance(visitInfo, lsst.afw.image.VisitInfo)
        # We should always get a valid date and exposure time
        self.assertIsInstance(visitInfo.getDate(), lsst.daf.base.DateTime)
        self.assertTrue(visitInfo.getDate().isValid())
        self.assertEqual(visitInfo.getExposureTime(), self.butler_get_data.exptimes["raw"])

    def testRawFilename(self):
        filepath = self.butler.get("raw_filename", self.dataIds["raw"])[0]
        if "[" in filepath:  # Remove trailing HDU specifier
            filepath = filepath[:filepath.find("[")]
        filename = os.path.split(filepath)[1]
        self.assertEqual(filename, self.mapper_data.raw_filename)

    def testQueryRawAmp(self):
        # Base the tests on the first reference metadata query
        formats = self.mapper_data.query_format.copy()
        query, expect = self.mapper_data.queryMetadata[0]
        result = self.mapper.queryMetadata("raw_amp", formats, query)
        self.assertEqual(sorted(result), sorted(expect))

        # Listing all channels -- we expect the result to be the expected
        # result copied for each channel
        formats = formats.copy()
        formats.insert(0, "channel")
        expectall = [(i+1, *expect[0]) for i in range(16)]
        result = self.mapper.queryMetadata("raw_amp", formats, query)
        self.assertEqual(sorted(result), sorted(expectall))

        # Now fix a channel
        query = query.copy()
        query["channel"] = 3
        result = self.mapper.queryMetadata("raw_amp", formats, query)
        expect = [(query["channel"], *expect[0])]
        self.assertEqual(sorted(result), sorted(expect))

        # Fix a channel out of range
        query["channel"] = 20
        with self.assertRaises(ValueError):
            self.mapper.queryMetadata("raw_amp", formats, query)

    def _testCoaddId(self, idName):
        coaddId = self.butler.get(idName, dataId={"tract": 9813, "patch": "3,4",
                                                  "filter": self.butler_get_data.filters["raw"]})
        self.assertIsInstance(coaddId, int)
        maxbits = self.butler.get(f"{idName}_bits")
        self.assertLess(maxbits, 64)
        self.assertLess(coaddId.bit_length(), maxbits, f"compare bit length for {idName}")

        # Check failure modes
        with self.assertRaises(RuntimeError):
            self.butler.get(idName, dataId={"tract": 9813, "patch": "-3,4"})
        with self.assertRaises(RuntimeError):
            self.butler.get(idName, dataId={"tract": -9813, "patch": "3,4"})
        with self.assertRaises(RuntimeError):
            self.butler.get(idName, dataId={"tract": 2**self.mapper._nbit_tract+1, "patch": "3,4"})
        with self.assertRaises(RuntimeError):
            self.butler.get(idName, dataId={"tract": 2, "patch": f"3,{2**self.mapper._nbit_patch+1}"})

    def testDeepCoaddId(self):
        self._testCoaddId("deepCoaddId")

    def testDcrCoaddId(self):
        self._testCoaddId("dcrCoaddId")

    def testDeepMergedCoaddId(self):
        self._testCoaddId("deepMergedCoaddId")

    def testDcrMergedCoaddId(self):
        self._testCoaddId("dcrMergedCoaddId")

    def testCcdExposureIdBits(self):
        """Check that we have enough bits for the exposure ID"""
        bits = self.butler.get('ccdExposureId_bits')
        ccdExposureId = self.butler_get_data.exposureIds["raw"]
        self.assertLessEqual(ccdExposureId.bit_length(), bits,
                             f"Can detector_exposure_id {ccdExposureId} fit in {bits} bits")
