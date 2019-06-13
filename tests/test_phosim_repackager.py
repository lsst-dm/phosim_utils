#
# LSST Data Management System
# Copyright 2012-2017 LSST Corporation.
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

import os
import shutil
import astropy.io.fits as fits
import unittest

from lsst.utils import getPackageDir
from lsst.phosim.utils.phosim_repackager import PhoSimRepackager


class TestPhoSimRepackager(unittest.TestCase):
    """Test the PhoSim Repackager class."""

    def setUp(self):

        self.phoSim_repackager = PhoSimRepackager()

        package_dir = getPackageDir("phosim_utils")
        test_dir = os.path.join(package_dir, "tests")

        self.tmp_test_dir = os.path.join(test_dir, "tmp")
        self._make_dir(self.tmp_test_dir)

        self.test_data_dir = os.path.join(test_dir, "testData")
        self.base_eimg_file_name = "lsst_e_9006001_f1_R22_S00_E000.fits"
        self.eimg_file_path = os.path.join(
            self.test_data_dir, "%s.gz" % self.base_eimg_file_name)

    def _make_dir(self, directory):

        if (not os.path.exists(directory)):
            os.makedirs(directory)

    def tearDown(self):

        shutil.rmtree(self.tmp_test_dir)

    def test_process_visit_eimage(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager.process_visit_eimage(
            self.test_data_dir, out_dir=self.tmp_test_dir, prefix='lsst')

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

    def test_repackage_eimage(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager.repackage_eimage(
            self.eimg_file_path, out_dir=self.tmp_test_dir)

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

        file_name = self.base_eimg_file_name
        self._check_repackaged_eimg_file(file_name)

    def _check_repackaged_eimg_file(self, file_name):

        eimg_file = os.path.join(self.tmp_test_dir, file_name)
        sensor = fits.open(eimg_file)[0]

        # Check the image is transposed
        self.assertEqual(sensor.data.shape, (4000, 4072))

        # Check the header imformation
        header = sensor.header
        self.assertEqual(header['RAFTNAME'], 'R22')
        self.assertEqual(header['SENSNAME'], 'S00')
        self.assertEqual(header['ROTANGLE'], 0.0)

    def _get_num_of_file_in_dir(self, directory):

        return len([name for name in os.listdir(directory)
                   if os.path.isfile(os.path.join(directory, name))])

    def test_mef_filename_eimage(self):

        out_dir = os.path.join(os.sep, "path", "of", "output")

        outfile = PhoSimRepackager.mef_filename_eimage(
            self.eimg_file_path, out_dir=out_dir)

        ans_outfile = os.path.join(out_dir, self.base_eimg_file_name)
        self.assertEqual(outfile, ans_outfile)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
