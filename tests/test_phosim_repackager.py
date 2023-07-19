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
from lsst.phosim.utils.phosim_repackager import (
    PhoSimRepackager,
    updateComCamSpecificData,
    updateLsstCamSpecificData
)


class TestPhoSimRepackager(unittest.TestCase):
    """Test the PhoSim Repackager class."""

    @classmethod
    def setUpClass(cls):
        cls.phoSim_repackager_lsstcam = PhoSimRepackager(
            instName="lsst", image_type="skyexp", focusz=-1500
        )
        cls.phoSim_repackager_calib = PhoSimRepackager(
            instName="lsst", image_type="dark"
        )
        cls.phoSim_repackager_comcam = PhoSimRepackager(
            "comcam", image_type="skyexp", focusz=-1500
        )
        cls.phoSim_repackager_comcam_derotate = PhoSimRepackager(
            "comcam", image_type="skyexp", focusz=-1500, derotate=False
        )

    def setUp(self):

        package_dir = getPackageDir("phosim_utils")
        test_dir = os.path.join(package_dir, "tests")

        self.tmp_test_dir = os.path.join(test_dir, "tmp")
        self._make_dir(self.tmp_test_dir)

        self.test_data_dir_eimg = os.path.join(test_dir, "testData", "eimage")
        self.test_data_dir_amp = os.path.join(test_dir, "testData", "amp")

        self.base_eimg_file_name = "lsst_e_9006001_f1_R22_S22_E000.fits"
        self.repackaged_eimg_file_name = "MC_H_20000217_006001_R22_S22_e.fits"
        self.repackaged_amp_file_name = "MC_H_20000217_006001_R22_S22.fits"

        self.base_eimg_file_name_comcam = "comcam_e_9006001_f1_R22_S20_E000.fits"
        self.repackaged_eimg_file_name_comcam = "CC_H_20000217_006001_R22_S20_e.fits"
        self.repackaged_amp_file_name_comcam = "CC_H_20000217_006001_R22_S20.fits"

        self.eimg_file_path = os.path.join(
            self.test_data_dir_eimg, "%s.gz" % self.base_eimg_file_name
        )
        self.amp_file_paths = sorted(
            glob.glob(os.path.join(self.test_data_dir_amp, "lsst_a_*"))
        )

        self.eimg_file_path_comcam = os.path.join(
            self.test_data_dir_eimg, "%s.gz" % self.base_eimg_file_name_comcam
        )
        self.amp_file_paths_comcam = sorted(
            glob.glob(os.path.join(self.test_data_dir_amp, "comcam_a_*"))
        )

    def _make_dir(self, directory):

        if not os.path.exists(directory):
            os.makedirs(directory)

    def tearDown(self):

        shutil.rmtree(self.tmp_test_dir)

    def test_phoSim_repackager_init(self):

        # controller is the same for both instruments
        self.assertEqual(self.phoSim_repackager_lsstcam.CONTRLLR, "H")
        self.assertEqual(self.phoSim_repackager_comcam.CONTRLLR, "H")
        self.assertEqual(self.phoSim_repackager_comcam_derotate.CONTRLLR, "H")

        # telescope code is different for each instrument
        self.assertEqual(self.phoSim_repackager_lsstcam.telcode, "MC")
        self.assertEqual(self.phoSim_repackager_comcam.telcode, "CC")
        self.assertEqual(self.phoSim_repackager_comcam_derotate.telcode, "CC")

        # instrument name is different for each instrument
        self.assertEqual(self.phoSim_repackager_lsstcam.instrument, "lsstCam")
        self.assertEqual(self.phoSim_repackager_comcam.instrument, "comCam")
        self.assertEqual(self.phoSim_repackager_comcam_derotate.instrument, "comCam")

        # image_type by default is skyexp
        self.assertEqual(self.phoSim_repackager_lsstcam.image_type, "skyexp")

        # image_type can be also dark, flat, bias
        self.assertEqual(self.phoSim_repackager_calib.image_type, "dark")

        # focusz is extra-focal, i.e -1500 for lsstcam and comcam test files
        # but 0 for calibs
        self.assertEqual(self.phoSim_repackager_lsstcam.focusz, -1500)
        self.assertEqual(self.phoSim_repackager_comcam.focusz, -1500)
        self.assertEqual(self.phoSim_repackager_comcam_derotate.focusz, -1500)
        self.assertEqual(self.phoSim_repackager_calib.focusz, 0)

        # Test derotate Flag
        self.assertTrue(self.phoSim_repackager_lsstcam.derotate)
        self.assertTrue(self.phoSim_repackager_comcam.derotate)
        self.assertFalse(self.phoSim_repackager_comcam_derotate.derotate)

    def test_process_visit_eimage(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager_lsstcam.process_visit_eimage(
            self.test_data_dir_eimg, out_dir=self.tmp_test_dir, instName="lsst"
        )

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

    def test_process_visit_eimage_comcam(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager_comcam.process_visit_eimage(
            self.test_data_dir_eimg, out_dir=self.tmp_test_dir, instName="comcam"
        )

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

    def test_repackage_eimage(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager_lsstcam.repackage_eimage(
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
        # and the dimension -
        # for lsstCam,  R22 is E2V
        self.assertEqual(sensor.data.shape, (4004, 4096))

        # Check the header information
        header = sensor.header
        self.assertEqual(header["RAFTBAY"], "R22")
        self.assertEqual(header["CCDSLOT"], "S22")
        self.assertEqual(header["RA"], 0.0)

        # Close the FITS file
        hdul.close()

    def test_repackage_eimage_comcam(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager_comcam.repackage_eimage(
            self.eimg_file_path_comcam, out_dir=self.tmp_test_dir
        )

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

        file_name = self.repackaged_eimg_file_name_comcam
        self._check_repackaged_eimg_file_comcam(file_name)

    def _check_repackaged_eimg_file_comcam(self, file_name):

        eimg_file = os.path.join(self.tmp_test_dir, file_name)
        hdul = fits.open(eimg_file)
        sensor = hdul[0]

        # Check the image is transposed
        # and the dimension correct -
        # for comcam, R22 is ITL
        self.assertEqual(sensor.data.shape, (4000, 4072))

        # Check the header imformation
        header = sensor.header
        self.assertEqual(header["RAFTBAY"], "R22")
        self.assertEqual(header["CCDSLOT"], "S20")
        self.assertEqual(header["RA"], 0.0)

        # Close the FITS file
        hdul.close()

    def test_process_visit_amp_image(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager_lsstcam.process_visit(
            self.test_data_dir_amp, out_dir=self.tmp_test_dir, instName="lsst"
        )

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

    def test_process_visit_amp_image_comcam(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager_comcam.process_visit(
            self.test_data_dir_amp, out_dir=self.tmp_test_dir, instName="comcam"
        )

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

    def test_repackage_amp_image(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager_lsstcam.repackage(
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
        # for lsstCam, R22 is E2V
        self.assertEqual(hdul[1].data.shape, (2048, 576))

        # Check the number of amps: main header
        # plus 16 amps
        self.assertEqual(len(hdul), 17)

        # Check the main header information
        header = hdul[0].header
        self.assertEqual(header["RAFTBAY"], "R22")
        self.assertEqual(header["CCDSLOT"], "S22")
        self.assertEqual(header["LSST_NUM"], "E2V-CCD250-370")
        self.assertEqual(header["OBSID"], "MC_H_20000217_006001")
        self.assertEqual(header["RA"], 0.0)
        self.assertEqual(header["FILTER"], "g_6")
        self.assertEqual(header["IMGTYPE"], "SKYEXP")
        self.assertEqual(header["FOCUSZ"], -1500)

        # Check the amp header information
        header = hdul[1].header
        self.assertEqual(header["AMPID"], "C10")
        self.assertEqual(header["CCDID"], "R22_S22")
        self.assertEqual(header["EXTNAME"], "Segment10")
        self.assertEqual(header["ROTANG"], 90)

        # Close the file
        hdul.close()

    def test_repackage_amp_image_comcam(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager_comcam.repackage(
            self.amp_file_paths_comcam, out_dir=self.tmp_test_dir
        )

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

        file_name = self.repackaged_amp_file_name_comcam
        self._check_repackaged_amp_file_comcam(
            file_name, derotate=self.phoSim_repackager_comcam.derotate
        )

    def test_repackage_amp_image_comcam_derotate(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager_comcam_derotate.repackage(
            self.amp_file_paths_comcam, out_dir=self.tmp_test_dir
        )

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

        file_name = self.repackaged_amp_file_name_comcam
        self._check_repackaged_amp_file_comcam(
            file_name, derotate=self.phoSim_repackager_comcam_derotate.derotate
        )

    def _check_repackaged_amp_file_comcam(self, file_name, derotate):

        amp_file = os.path.join(self.tmp_test_dir, file_name)
        hdul = fits.open(amp_file)

        # Check the amp has correct shape
        # for comcam, R22 is ITL
        self.assertEqual(hdul[1].data.shape, (2048, 544))

        # Check the number of amps: main header
        # plus 16 amps
        self.assertEqual(len(hdul), 17)

        # Check the main header information
        header = hdul[0].header
        self.assertEqual(header["RAFTBAY"], "R22")
        self.assertEqual(header["CCDSLOT"], "S20")
        self.assertEqual(header["LSST_NUM"], "ITL-3800C-319")
        self.assertEqual(header["OBSID"], "CC_H_20000217_006001")
        self.assertEqual(header["RA"], 0.0)
        self.assertEqual(header["FILTER"], "g_01")
        self.assertEqual(header["IMGTYPE"], "SKYEXP")
        self.assertEqual(header["FOCUSZ"], -1500)

        # Check the amp header information
        header = hdul[1].header
        self.assertEqual(header["AMPID"], "C10")
        self.assertEqual(header["CCDID"], "R22_S20")
        self.assertEqual(header["EXTNAME"], "Segment10")
        self.assertEqual(header["ROTANG"], 90 if derotate else 0)

        # Close the file
        hdul.close()

    def test_repackage_amp_calib_lsst(self):

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 0)

        self.phoSim_repackager_calib.repackage(
            self.amp_file_paths, out_dir=self.tmp_test_dir
        )

        num_file = self._get_num_of_file_in_dir(self.tmp_test_dir)
        self.assertEqual(num_file, 1)

        file_name = self.repackaged_amp_file_name
        self._check_repackaged_amp_calib_lsst(file_name)

    def _check_repackaged_amp_calib_lsst(self, file_name):

        amp_file = os.path.join(self.tmp_test_dir, file_name)
        hdul = fits.open(amp_file)

        # Check the amp has correct shape
        # for lsst, R22 is ITL
        self.assertEqual(hdul[1].data.shape, (2048, 576))

        # Check the number of amps: main header
        # plus 16 amps
        self.assertEqual(len(hdul), 17)

        # Check the main header information
        header = hdul[0].header
        self.assertEqual(header["RAFTBAY"], "R22")
        self.assertEqual(header["CCDSLOT"], "S22")
        self.assertEqual(header["LSST_NUM"], "E2V-CCD250-370")
        self.assertEqual(header["OBSID"], "MC_H_20000217_006001")
        self.assertEqual(header["RA"], 0.0)
        self.assertEqual(header["FILTER"], "g_6")
        self.assertEqual(header["IMGTYPE"], "DARK")
        self.assertEqual(header["FOCUSZ"], 0)

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

    def test_updateComCamSpecificData(self):

        # update header for raw comcam amp file
        hdul = fits.open(self.amp_file_paths_comcam[0])
        input_header = hdul[0].header
        self.assertEqual(input_header["FILTER"], "g")

        # test that the header got updated
        updated_header = updateComCamSpecificData(input_header)
        self.assertEqual(updated_header["FILTER"], "g_01")

        # Close the file
        hdul.close()

        # try updating header with unknown filter name
        wrong_header = fits.Header({"FILTER": "x"})
        self.assertRaises(ValueError, updateComCamSpecificData, wrong_header)

    def test_updateLsstCamSpecificData(self):

        # update header for raw comcam amp file
        hdul = fits.open(self.amp_file_paths[0])
        input_header = hdul[0].header
        self.assertEqual(input_header["FILTER"], "g")

        # test that the header got updated
        updated_header = updateLsstCamSpecificData(input_header)
        self.assertEqual(updated_header["FILTER"], "g_6")

        # Close the file
        hdul.close()

        # try updating header with unknown filter name
        wrong_header = fits.Header({"FILTER": "x"})
        self.assertRaises(ValueError, updateLsstCamSpecificData, wrong_header)


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
