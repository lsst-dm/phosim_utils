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
import glob

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

        self.test_data_dir_eimg = os.path.join(test_dir, "testData", "eimage")
        self.test_data_dir_amp = os.path.join(test_dir, "testData", "amp")
        self.base_eimg_file_name = "lsst_e_9006001_f1_R22_S22_E000.fits"
        self.repackaged_eimg_file_name = "MC_H_20000217_006001_R22_S22_e.fits"
        self.repackaged_amp_file_name = "MC_H_20000217_006001_R22_S22.fits"
        self.eimg_file_path = os.path.join(
            self.test_data_dir_eimg, "%s.gz" % self.base_eimg_file_name)
        self.amp_file_paths = sorted(glob.glob(os.path.join(self.test_data_dir_amp,
                                     'lsst_a_*')))

    def _make_dir(self, directory):

        if not os.path.exists(directory):
            os.makedirs(directory)

    def tearDown(self):

        shutil.rmtree(self.tmp_test_dir)

    def test_process_visit_eimage(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager.process_visit_eimage(
            self.test_data_dir_eimg, out_dir=self.tmp_test_dir, instName="lsst"
        )

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

    def test_repackage_eimage(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager.repackage_eimage(
            self.eimg_file_path, out_dir=self.tmp_test_dir
        )

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

        file_name = self.repackaged_eimg_file_name
        self._check_repackaged_eimg_file(file_name)

    def _check_repackaged_eimg_file(self, file_name):

        eimg_file = os.path.join(self.tmp_test_dir, file_name)
        hdul = fits.open(eimg_file)
        sensor = hdul[0]

        # Check the image is transposed
        self.assertEqual(sensor.data.shape, (4004, 4096))

        # Check the header imformation
        header = sensor.header
        self.assertEqual(header["RAFTBAY"], "R22")
        self.assertEqual(header["CCDSLOT"], "S22")
        self.assertEqual(header["RA"], 0.0)

        # Close the FITS file
        hdul.close()

    def test_process_visit_amp_image(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager.process_visit(
            self.test_data_dir_amp, out_dir=self.tmp_test_dir, instName="lsst"
        )

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

    def test_repackage_amp_image(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager.repackage(
            self.amp_file_paths, out_dir=self.tmp_test_dir
        )

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

        file_name = self.repackaged_amp_file_name
        self._check_repackaged_amp_file(file_name)

    def _check_repackaged_amp_file(self, file_name):

        amp_file = os.path.join(self.tmp_test_dir, file_name)
        hdul = fits.open(amp_file)

        # Check the amp has correct shape
        self.assertEqual(hdul[1].data.shape, (2048, 576))

        # Check the number of amps: main header
        # plus 16 amps
        self.assertEqual(len(hdul), 17)

        # Check the main header imformation
        header = hdul[0].header
        self.assertEqual(header["RAFTBAY"], "R22")
        self.assertEqual(header["CCDSLOT"], "S22")
        self.assertEqual(header["LSST_NUM"], "E2V-CCD250-370")
        self.assertEqual(header["OBSID"], "MC_H_20000217_006001")
        self.assertEqual(header["RA"], 0.0)

        # Check the amp header information
        header = hdul[1].header
        self.assertEqual(header["AMPID"], "C10")
        self.assertEqual(header["CCDID"], "R22_S22")
        self.assertEqual(header["EXTNAME"], "Segment10")

        # Close the file
        hdul.close()

    def _get_num_of_file_in_dir(self, directory):

        return len(
            [
                name
                for name in os.listdir(directory)
                if os.path.isfile(os.path.join(directory, name))
            ]
        )


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
