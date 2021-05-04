"""
Code to repackage phosim files. The amplifier files will be transformed into
single sensor multi-extension FITS files, with an HDU per amplifier. For the
eimage, the header will be updated to match the translator's definition in
obs_lsst.
"""
import os
import sys
import glob
import time
from collections import defaultdict
import astropy.io.fits as fits
import astropy.time
from lsst.obs.lsst import LsstCam, LsstComCam

__all__ = ["PhoSimRepackager", "noao_section_keyword"]


def noao_section_keyword(bbox, flipx=False, flipy=False):
    """
    Convert bounding boxes into NOAO section keywords.

    Parameters
    ----------
    bbox : lsst.afw.geom.Box2I
        Bounding box.
    flipx : bool, optional
        Flag to indicate that data should be flipped in the x-direction.
        (default is False)
    flipy : bool, optional
        Flag to indicate that data should be flipped in the y-direction.
        (default is False)
    """
    xmin, xmax = bbox.getMinX() + 1, bbox.getMaxX() + 1
    ymin, ymax = bbox.getMinY() + 1, bbox.getMaxY() + 1
    if flipx:
        xmin, xmax = xmax, xmin
    if flipy:
        ymin, ymax = ymax, ymin
    return "[%i:%i,%i:%i]" % (xmin, xmax, ymin, ymax)


def updateComCamSpecificData(header):
    """
    Update comCam specific data.
    Currently translate phosim
    filter name (one of u,g,r,i,z,y) to that appropriate
    for obs_lsst - see DM-21706.

    Parameters
    ----------
    header : FITS header of a phosim output files.

    Returns
    --------
    header : updated FITS header.

    Raises
    ------
    ValueError
        Not possible to update the header due to
        incorrect filter name.
    """
    phosim_to_comcam = {
        "u": "06",
        "g": "01",
        "r": "03",
        "i": "06",
        "z": "02",
        "y": "04",
    }

    phosimFilter = header["FILTER"]

    # test if filter can be translated
    if phosimFilter in phosim_to_comcam.keys():
        suffix = phosim_to_comcam[phosimFilter]
        header["FILTER"] = f"{phosimFilter}_{suffix}"
        return header  # return updated FITS header
    else:
        raise ValueError(
            f"This phosim filter:{phosimFilter} cannot \
            be translated for comcam")


class PhoSimRepackager:
    """
    Class to repackage phosim amplifier files into single sensor
    MEFs with one HDU per amp, using LSE-400 for header conformity.
    """

    CONTRLLR = "H"  # Controller: H stands for pHosim, to distinguish it from the
    # observed (not simulated) data,  added via DM-27863
    # to obs_lsst translators/lsst.py

    def __init__(self, instName="lsst", image_type="skyexp", focusz=0):
        """
        Parameters
        ----------
        instName: str, optional
            Instrument to use: comcam or lsst. Corresponds
            to the prefix of phosim amplifier files.
            Note: corner sensors are included in lsst.
            (the default is 'lsst').
        image_type: str, optional
            Image type  corresponds to the IMGTYPE header,
            those approved by obs_lsst metadata translator
            for lsstCam or comCam are SKYEXP, FLAT, DARK,
            BIAS. (the default is 'skyexp').
        focusz: int, optional
            The position of the main camera hexapod in micrometers.
            For in-focus image it is 0.
            Assuming defocal distance of d micrometers,
            for extra-focal image it would be -d,
            and for intra-focal image it would be d.
            Added to the repackaged image
            primary image header as FOCUSZ.
            (the default is 0).

        (the default is 0)
        """
        # Use appropriate obs_lsst mapper camera object
        # and telescope code, and declare the image type
        # to be stored in the  header.
        self.image_type = image_type
        self.focusz = focusz
        if instName == "lsst":
            self.camera = LsstCam().getCamera()
            self.telcode = "MC"  # Main Camera
            self.instrument = "lsstCam"

        elif instName == "comcam":
            self.camera = LsstComCam().getCamera()
            self.telcode = "CC"  # Commissioning Camera
            self.instrument = "comCam"
        else:
            raise ValueError(f"This instrument name ({instName}) is not supported.")

    def process_visit(self, visit_dir, out_dir=None, instName="lsst", verbose=False):
        """
        Parameters
        ----------
        visit_dir: str
            Directory containing the phosim amplifier for a given visit.
        out_dir: str, optional
            Output directory for MEF files. If None, then a directory
            with name v<visit #>-<band> will be created in the cwd.
            (the default is None)
        instName: str, optional
            Name of instrument to use: either comcam or lsst.
            Corresponds to the prefix of phosim amplifier files in visit_dir.
            (the default is 'lsst')
        verbose: bool, optional
            Set to True to print out time for processing each sensor.
            (the default is False)
        """
        phosim_amp_files = sorted(glob.glob(os.path.join(visit_dir, f"{instName}_a_*")))
        amp_files = defaultdict(list)

        for item in phosim_amp_files:
            sensor_id = "_".join(os.path.basename(item).split("_")[4:6])
            amp_files[sensor_id].append(item)

        if out_dir is None:
            tokens = os.path.basename(phosim_amp_files[0]).split("_")
            out_dir = "v%07i-%s" % (int(tokens[2]), "ugrizy"[int(tokens[3][1])])

        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        for sensor_id in amp_files:
            # print time spent repackaging
            t0 = time.time()

            # repackage with the same function for corner or main raft sensors
            self.repackage(amp_files[sensor_id], out_dir=out_dir)

            if verbose:
                ttl_time = time.time() - t0
                print(f"Repackaging took {ttl_time:.3f} seconds")
                sys.stdout.flush()

    def repackage(self, phosim_amp_files, out_dir="."):
        """
        Repackage a collection of phosim amplifier files for a
        single sensor into a multi-extension FITS file.

        Parameters
        ----------
        phosim_amp_files: list
            List of phosim amplifier filenames.
        out_dir: str, optional
            The path to output directory. (the default is '.')
        """

        segments = {}
        for filename in phosim_amp_files:
            phosim_channel = os.path.basename(filename).split("_")[6][1:]
            segments[phosim_channel] = fits.open(filename)[0]

        # Read raft name from the first filename
        raft = os.path.basename(filename).split("_")[4]
        ccdslot = os.path.basename(filename).split("_")[5]
        sensorId = "%s_%s" % (raft, ccdslot)
        print(f"Repackaging sensor {sensorId}")

        detectors = self.camera.get(sensorId)

        # initialize the FITS file
        sensor = fits.HDUList(fits.PrimaryHDU())

        # Set the NOAO section keywords based on the pixel geometry
        # in the obs_lsst object for each amplifier header
        for amp in detectors:
            ampName = amp.getName()[1:]
            hdu = segments[ampName]
            hdu.header["EXTNAME"] = "Segment%s" % ampName
            hdu.header["DATASEC"] = noao_section_keyword(amp.getRawDataBBox())
            hdu.header["DETSEC"] = noao_section_keyword(
                amp.getBBox(), flipx=amp.getRawFlipX(), flipy=amp.getRawFlipY()
            )
            # Offset the rotation angle by 90 degrees
            # to set correct WCS
            hdu.header["ROTANG"] = 90 - hdu.header["ROTANG"]

            # Set the filter information.
            # For lsstCam, the phosim filter names (ugrizy)
            # are appropriate, but for comCam they need
            # to have a different name as in obs_lsst filters.py
            if self.instrument == "comCam":
                hdu.header = updateComCamSpecificData(hdu.header)

            # Remove the incorrect BIASSEC keyword that phosim writes
            try:
                hdu.header.remove("BIASSEC")
            except KeyError:
                pass
            sensor.append(hdu)

        # Set keywords in primary HDU, extracting most of the relevant
        # ones from the first phosim amplifier file.
        sensor[0].header["EXPTIME"] = sensor[1].header["EXPTIME"]
        sensor[0].header["DARKTIME"] = sensor[1].header["DARKTIME"]

        # Add keyword for defocal position, FOCUSZ
        sensor[0].header["FOCUSZ"] = self.focusz
        # Call phosim OBSID as RUNNUM to be handled properly
        # by obs_lsst LsstCam or LsstComCam translators
        sensor[0].header["RUNNUM"] = sensor[1].header["OBSID"]
        sensor[0].header["MJD-OBS"] = sensor[1].header["MJD-OBS"]
        DATEOBS = astropy.time.Time(sensor[1].header["MJD-OBS"], format="mjd").isot
        sensor[0].header["DATE-OBS"] = DATEOBS
        YEAR, MONTH, DAYTIME = DATEOBS.split("-")
        DAY = DAYTIME[:2]
        DAYOBS = "%s%s%s" % (YEAR, MONTH, DAY)
        sensor[0].header["DAYOBS"] = DAYOBS

        DATE = sensor[1].header["DATE"]  # file creation date
        sensor[0].header["DATE"] = DATE
        sensor[0].header["MJD"] = astropy.time.Time(DATE, format="isot").mjd

        sensor[0].header["FILTER"] = sensor[1].header["FILTER"]

        serial = detectors.getSerial()  # eg. ITL-4400B-029
        CCD_MANU, CCD_TYPE, CCD_NUM = serial.split("-")
        sensor[0].header["LSST_NUM"] = serial
        sensor[0].header["CCD_MANU"] = CCD_MANU  # eg. ITL
        sensor[0].header["CCD_TYPE"] = CCD_TYPE  # eg. 4400B

        sensor[0].header["DETSIZE"] = noao_section_keyword(detectors.getBBox())
        sensor[0].header["INSTRUME"] = self.instrument
        sensor[0].header["TELESCOP"] = "LSST"

        sensor[0].header["TELCODE"] = self.telcode
        sensor[0].header["CONTRLLR"] = self.CONTRLLR

        # Set the integer sequence number. Since it can be 6 digits long,
        # construct from phosim OBSID - eg. 9006001, making an integer from
        # the last 6 digits of OBSID,  i.e. 6001
        SEQNUM = int(sensor[0].header["RUNNUM"][-6:])
        sensor[0].header["SEQNUM"] = SEQNUM

        # Contruct new OBSID from telescope code (eg. MC),
        # controller code (eg. H), observation date (eg. 20210325),
        # sequence number (eg. 6001), which is zero-padded
        # to match what happens in real CCS
        OBSID = f"{self.telcode}_{self.CONTRLLR}_{DAYOBS}_{SEQNUM:06d}"
        sensor[0].header["OBSID"] = OBSID

        sensor[0].header["TESTTYPE"] = "PHOSIM"
        sensor[0].header["IMGTYPE"] = self.image_type.upper()  # "SKYEXP", "BIAS", "FLAT", "DARK"
        sensor[0].header["RAFTBAY"] = raft
        sensor[0].header["CCDSLOT"] = ccdslot
        sensor[0].header["RASTART"] = sensor[1].header["RA_DEG"]
        sensor[0].header["DECSTART"] = sensor[1].header["DEC_DEG"]
        sensor[0].header["ROTSTART"] = sensor[1].header["ROTANG"]
        sensor[0].header["RA"] = sensor[1].header["RA_DEG"]
        sensor[0].header["DEC"] = sensor[1].header["DEC_DEG"]
        sensor[0].header["ROTPA"] = sensor[1].header["ROTANG"]
        sensor[0].header["ROTCOORD"] = "sky"
        sensor[0].header["AMSTART"] = sensor[1].header["AIRMASS"]
        sensor[0].header["ELSTART"] = (
            90.0 - sensor[1].header["ZENITH"]
        )  # altitude angle

        # Filename created from TELCODE, CONTRLLR, DAYOBS, SEQNUM , raft, ccdslot.
        # eg. MC_C_20200825_000032_R00_SW0.fits
        # Below, OBSID already contains TELCODE, CONTRLLR, DAYOBS, str(SEQNUM).zfill(6)
        filename = "%s_%s_%s.fits" % (OBSID, raft, ccdslot)
        filename = os.path.join(out_dir, filename)

        # save the FITS file
        sensor.writeto(filename, overwrite=True)
        print(f"Saved as {filename}")

    def process_visit_eimage(
        self, visit_dir, out_dir=None, instName="lsst", verbose=False
    ):
        """Process the visit eimage data.

        Parameters
        ----------
        visit_dir: str
            Directory containing the phosim eimage for a given visit.
        out_dir : str, optional
            Output directory for repackaged files. If None, then a directory
            with name 'eimage' will be created in the cwd.
            (the default is None.)
        instName: str, optional
            Name of instrument to use: either comcam or lsst.
            Corresponds to the prefix of phosim eimage files in visit_dir.
            (the default is 'lsst')
        verbose : bool, optional
            Set to True to print out time for processing each sensor.
            (the default is False.)
        """

        phosim_eimg_files = sorted(
            glob.glob(os.path.join(visit_dir, f"{instName}_e_*"))
        )

        if out_dir is None:
            out_dir = "eimage"

        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        for eimg_file in phosim_eimg_files:
            if verbose:
                sys.stdout.write(eimg_file + "  ")
            t0 = time.time()

            print(f"Repackaging {eimg_file}")
            self.repackage_eimage(eimg_file, out_dir=out_dir)
            if verbose:
                print(time.time() - t0)
                sys.stdout.flush()

    def repackage_eimage(self, phosim_eimg_file, out_dir="."):
        """Repackage the phosim eimage file to the format required by obs_lsst.

        Parameters
        ----------
        phosim_eimg_file : str
            PhoSim eimage file.
        out_dir : str, optional
            Output directory (the default is '.')
        """

        sensor = fits.open(phosim_eimg_file)[0]

        raft = os.path.basename(phosim_eimg_file).split("_")[4]
        ccdslot = os.path.basename(phosim_eimg_file).split("_")[5]
        sensorId = "%s_%s" % (raft, ccdslot)
        print(f"Repackaging sensor {sensorId}")

        detectors = self.camera.get(sensorId)

        # Set BOT-like keywords in the HDU
        sensor.header["RUNNUM"] = sensor.header["OBSID"]

        DATEOBS = astropy.time.Time(sensor.header["MJD-OBS"], format="mjd").isot
        sensor.header["DATE-OBS"] = DATEOBS

        YEAR, MONTH, DAYTIME = DATEOBS.split("-")
        DAY = DAYTIME[:2]
        DAYOBS = "%s%s%s" % (YEAR, MONTH, DAY)
        sensor.header["DAYOBS"] = DAYOBS

        # file creation date as mjd
        sensor.header["MJD"] = astropy.time.Time(
            sensor.header["DATE"], format="isot"
        ).mjd

        # Set the filter information.
        # For lsstCam, the phosim filter names (ugrizy)
        # are appropriate, but for comCam they need
        # to have a different name as in obs_lsst filters.py
        if self.instrument == "comCam":
            sensor.header = updateComCamSpecificData(sensor.header)

        serial = detectors.getSerial()  # eg. ITL-4400B-029
        CCD_MANU, CCD_TYPE, CCD_NUM = serial.split("-")
        sensor.header["LSST_NUM"] = serial
        sensor.header["CCD_MANU"] = CCD_MANU  # eg. ITL
        sensor.header["CCD_TYPE"] = CCD_TYPE  # eg. 4400B

        sensor.header["DETSIZE"] = noao_section_keyword(detectors.getBBox())
        sensor.header["INSTRUME"] = self.instrument
        sensor.header["TELESCOP"] = "LSST"

        sensor.header["TELCODE"] = self.telcode
        sensor.header["CONTRLLR"] = self.CONTRLLR

        # Set the integer sequence number. Since it can be 6 digits long,
        # construct from phosim OBSID - eg. 9006001, making an integer from
        # the last 6 digits of OBSID,  i.e. 6001
        SEQNUM = int(sensor.header["RUNNUM"][-6:])
        sensor.header["SEQNUM"] = SEQNUM

        # Contruct new OBSID from telescope code (eg. MC),
        # controller code (eg. H), observation date (eg. 20210325),
        # sequence number (eg. 6001), which is zero-padded
        # to match what happens in real CCS
        OBSID = f"{self.telcode}_{self.CONTRLLR}_{DAYOBS}_{SEQNUM:06d}"
        sensor.header["OBSID"] = OBSID

        # offset by 90 degrees to match WCS
        sensor.header["ROTANG"] = 90 - sensor.header["ROTANG"]
        sensor.header["TESTTYPE"] = "PHOSIM"
        sensor.header["IMGTYPE"] = "SKYEXP"

        sensor.header["RAFTBAY"] = raft
        sensor.header["CCDSLOT"] = ccdslot

        sensor.header["RASTART"] = sensor.header["RA_DEG"]
        sensor.header["DECSTART"] = sensor.header["DEC_DEG"]
        sensor.header["ROTSTART"] = sensor.header["ROTANG"]

        sensor.header["RA"] = sensor.header["RA_DEG"]
        sensor.header["DEC"] = sensor.header["DEC_DEG"]
        sensor.header["ROTPA"] = sensor.header["ROTANG"]

        sensor.header["ROTANGLE"] = sensor.header["ROTANG"]
        sensor.header["ROTCOORD"] = "sky"

        sensor.header["AMSTART"] = sensor.header["AIRMASS"]
        sensor.header["ELSTART"] = 90.0 - sensor.header["ZENITH"]

        # Transpose the image to fulfill the geometry of postISRCCD by obs_lsst
        sensor.data = sensor.data.T

        # get a filename from TELCODE, CONTRLLR,  DAYOBS, SEQNUM , raft, ccdslot
        # eg. MC_H_20200825_000032_R00_SW0.fits
        filename = "%s_%s_%s_e.fits" % (OBSID, raft, ccdslot)
        filename = os.path.join(out_dir, filename)

        # save the FITS file
        sensor.writeto(filename, overwrite=True)
        print(f"Saved as {filename}")
