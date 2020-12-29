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
import numpy as np
from collections import defaultdict
import astropy.io.fits as fits
import astropy.time
from lsst.obs.lsst import LsstCam

__all__ = ['PhoSimRepackager', 'noao_section_keyword']


def noao_section_keyword(bbox, flipx=False, flipy=False):
    """
    Convert bounding boxes into NOAO section keywords.

    Parameters
    ----------
    bbox : lsst.afw.geom.Box2I
        Bounding box.
    flipx : bool
        Flag to indicate that data should be flipped in the x-direction.
    flipy : bool
        Flag to indicate that data should be flipped in the y-direction.
    """
    xmin, xmax = bbox.getMinX()+1, bbox.getMaxX()+1
    ymin, ymax = bbox.getMinY()+1, bbox.getMaxY()+1
    if flipx:
        xmin, xmax = xmax, xmin
    if flipy:
        ymin, ymax = ymax, ymin
    return '[%i:%i,%i:%i]' % (xmin, xmax, ymin, ymax)


class PhoSimRepackager:
    """
    Class to repackage phosim amplifier files into single sensor
    MEFs with one HDU per amp, using LSE-400 for header conformity.
    """
    def __init__(self):
        self.telcode = 'MC'
        self.contrllr = 'H'
        self.camera = LsstCam().getCamera()

    def process_visit(self, visit_dir, out_dir=None, prefix='lsst', verbose=False):
        """
        Parameters
        ----------
        visit_dir: str
            Directory containing the phosim amplifier for a given visit.
        out_dir: str [None]
            Output directory for MEF files. If None, then a directory
            with name v<visit #>-<band> will be created in the cwd.
        prefix: str ['lsst']
            Prefix to use to find PhoSim amplifier files in visit_dir.
        verbose: bool [False]
            Set to True to print out time for processing each sensor.
        """
        phosim_amp_files \
            = sorted(glob.glob(os.path.join(visit_dir, f'{prefix}_a_*')))
        amp_files = defaultdict(list)

        for item in phosim_amp_files:
            sensor_id = '_'.join(os.path.basename(item).split('_')[4:6])
            amp_files[sensor_id].append(item)

        if out_dir is None:
            tokens = os.path.basename(phosim_amp_files[0]).split('_')
            out_dir = 'v%07i-%s' % (int(tokens[2]), 'ugrizy'[int(tokens[3][1])])

        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        for sensor_id in amp_files:
            # print time spent repackaging
            if verbose:
                sys.stdout.write(sensor_id + '  ')
            t0 = time.time()

            # decide whether corner sensors or main raft
            if sensor_id in ['R00_S22', 'R04_S20', 'R40_S02', 'R44_S00']:
                print("repackaging corner")
                self.repackage_corner(amp_files[sensor_id], out_dir=out_dir,
                                      verbose=verbose)

            else:
                print("repackaging main")
                self.repackage_main(amp_files[sensor_id], out_dir=out_dir,
                                    verbose=verbose)

            if verbose:
                print(time.time() - t0)
                sys.stdout.flush()

    def repackage_corner(self, phosim_amp_files, out_dir='.', verbose=False):
        """
        Repackage a collection of phosim amplifier files for a
        single wavefront sensor into a two multi-extension FITS files,
        one for intra- and one for extra-focal segments.

        Parameters
        ----------
        phosim_amp_files: list
            List of phosim amplifier filenames.
        camera: lsst.afw.cameraGeom.camera.camera.Camera object
        out_dir : str, optional
            Output directory (the default is '.'.)
        verbose: bool, optional
            Set to True to print out processing information.
        """
        # PhoSim labels the amplifier channels differently than
        # LsstCamMapper . This is a mapping from the phosim channels
        # to the LsstCam channels.

        # SW1, intra-focal
        ch_map_intra = {'00': '17',
                        '01': '16',
                        '02': '15',
                        '03': '14',
                        '04': '13',
                        '05': '12',
                        '06': '11',
                        '07': '10'}

        # SW0, extra-focal
        ch_map_extra = {'10': '17',
                        '11': '16',
                        '12': '15',
                        '13': '14',
                        '14': '13',
                        '15': '12',
                        '16': '11',
                        '17': '10'}

        # Read the intra- and extra-focal amps as segments that will be
        # assembled into a fits file
        segmentsIntra = {}
        segmentsExtra = {}

        # Store intra- and extra-focal amps as keys of a segments
        # dictionary. Find channel name from phosim filename.
        # Correct the channel name into obs_lsst convention,
        # handling separately intra- and extra-focal cases.
        for filename in phosim_amp_files:
            phosim_channel = os.path.basename(filename).split('_')[6][1:]

            if phosim_channel in ch_map_intra.keys():
                segmentsIntra[ch_map_intra[phosim_channel]] = fits.open(filename)[0]

            elif phosim_channel in ch_map_extra.keys():
                segmentsExtra[ch_map_extra[phosim_channel]] = fits.open(filename)[0]

        # Read raft name from the last filename
        raft = os.path.basename(filename).split('_')[4]
        print('Repackaging sensor %s'%raft)

        TELCODE = self.telcode
        CONTRLLR = self.contrllr

        for ccdslot, segments in zip(['SW1', 'SW0'],
                                     [segmentsIntra, segmentsExtra]):
            print('\n', ccdslot, segments.keys())

            # initialize the FITS file
            sensor = fits.HDUList(fits.PrimaryHDU())
            sensorId = '%s_%s'%(raft, ccdslot)
            if verbose:
                print('Using amp info for %s'%sensorId)
            detectors = self.camera.get(sensorId)

            # Set the NOAO section keywords based on the pixel geometry
            # in the obs_lsst object for each amplifier header
            for amp in detectors:
                ampName = amp.getName()[1:]
                if verbose:
                    print(ampName)
                hdu = segments[ampName]
                hdu.header['EXTNAME'] = 'Segment%s' % ampName
                hdu.header['DATASEC'] = noao_section_keyword(amp.getRawDataBBox())
                hdu.header['DETSEC'] = noao_section_keyword(
                    amp.getBBox(),
                    flipx=amp.getRawFlipX(),
                    flipy=amp.getRawFlipY())

                # flip up-down the intra-focal images for correct orientation
                if ccdslot == 'SW1':
                    if verbose:
                        print('Flipping up-down..')
                    data = np.copy(hdu.data)
                    hdu.data[:] = np.flipud(data)

                # Remove the incorrect BIASSEC keyword that phosim writes
                try:
                    hdu.header.remove('BIASSEC')
                except KeyError:
                    pass
                sensor.append(hdu)

            # Set keywords in primary HDU, extracting most of the relevant
            # ones from the first phosim amplifier file.
            sensor[0].header['EXPTIME'] = sensor[1].header['EXPTIME']
            sensor[0].header['DARKTIME'] = sensor[1].header['DARKTIME']
            sensor[0].header['RUNNUM'] = sensor[1].header['OBSID']
            sensor[0].header['MJD-OBS'] = sensor[1].header['MJD-OBS']
            DATEOBS = astropy.time.Time(sensor[1].header['MJD-OBS'], format='mjd').isot
            sensor[0].header['DATE-OBS'] = DATEOBS
            YEAR, MONTH, DAYTIME = DATEOBS.split('-')
            DAY = DAYTIME[:2]
            DAYOBS = '%s%s%s'%(YEAR, MONTH, DAY)
            sensor[0].header['DAYOBS'] = DAYOBS

            DATE = sensor[1].header['DATE']  # file creation date
            sensor[0].header['DATE'] = DATE
            sensor[0].header['MJD'] = astropy.time.Time(DATE, format='isot').mjd

            sensor[0].header['FILTER'] = sensor[1].header['FILTER']

            serial = detectors.getSerial()  # eg. ITL-4400B-029
            CCD_MANU, CCD_TYPE, CCD_NUM = serial.split('-')
            sensor[0].header['LSST_NUM'] = serial
            sensor[0].header['CCD_MANU'] = CCD_MANU  # eg. ITL
            sensor[0].header['CCD_TYPE'] = CCD_TYPE  # eg. 4400B

            # NB: I get [1:4072,1:2000], whereas  BOT has '[1:4072,1:4000]'
            sensor[0].header['DETSIZE'] = noao_section_keyword(detectors.getBBox())
            sensor[0].header['INSTRUME'] = 'lsstCam'
            sensor[0].header['TELESCOP'] = 'LSST'

            sensor[0].header['TELCODE'] = TELCODE
            sensor[0].header['CONTRLLR'] = CONTRLLR
            # Set sequence number from OBSID, eg. 9006001, taking
            # the last 6 digits
            SEQNUM = int(sensor[0].header['OBSID'][-6:])
            sensor[0].header['SEQNUM'] = SEQNUM

            OBSID = "%s_%s_%s_%s"%(TELCODE, CONTRLLR, DAYOBS, str(SEQNUM).zfill(6))
            sensor[0].header['OBSID'] = OBSID

            sensor[0].header['TESTTYPE'] = 'PHOSIM'
            sensor[0].header['IMGTYPE'] = 'SKYEXP'
            sensor[0].header['RAFTBAY'] = raft
            sensor[0].header['CCDSLOT'] = ccdslot
            sensor[0].header['RASTART'] = sensor[1].header['RA_DEG']
            sensor[0].header['DECSTART'] = sensor[1].header['DEC_DEG']
            sensor[0].header['ROTSTART'] = sensor[1].header['ROTANG']
            sensor[0].header['RA'] = sensor[1].header['RA_DEG']
            sensor[0].header['DEC'] = sensor[1].header['DEC_DEG']
            sensor[0].header['ROTPA'] = sensor[1].header['ROTANG']
            sensor[0].header['ROTPOS'] = sensor[1].header['ROTANG']

            # Filename created from TELCODE, CONTRLLR, DAYOBS, SEQNUM , raft, ccdslot.
            # eg. MC_C_20200825_000032_R00_SW0.fits
            # Below, OBSID already contains TELCODE, CONTRLLR, DAYOBS, str(SEQNUM).zfill(6)
            filename = '%s_%s_%s.fits'%(OBSID, raft, ccdslot)
            filename = os.path.join(out_dir, filename)

            # save the FITS file
            sensor.writeto(filename, overwrite=True)
            print('Saved as %s'%filename)

    def repackage_main(self, phosim_amp_files, out_dir='.', verbose=False):
        """
        Repackage a collection of phosim amplifier files for a
        single wavefront sensor into a two multi-extension FITS files,
        one for intra- and one for extra-focal segments.

        Parameters
        ----------
        phosim_amp_files: list
            List of phosim amplifier filenames.
        out_dir : str, optional
            Output directory (the default is '.'.)
        verbose: bool, optional
            Set to True to print out processing information.
        """
        # Read the amplifiers into segments that will be
        # assembled into a fits file
        segments = {}
        for filename in phosim_amp_files:
            phosim_channel = os.path.basename(filename).split('_')[6][1:]
            segments[phosim_channel] = fits.open(filename)[0]

        # Read raft name from the last filename
        raft = os.path.basename(filename).split('_')[4]
        ccdslot = os.path.basename(filename).split('_')[5]
        print('Repackaging sensor %s'%raft)

        TELCODE = self.telcode
        CONTRLLR = self.contrllr

        ####
        # make  the FITS file
        ###
        if verbose:
            print('\n', ccdslot, segments.keys())

        # initialize the FITS file
        sensor = fits.HDUList(fits.PrimaryHDU())
        sensorId = '%s_%s'%(raft, ccdslot)
        if verbose:
            print('Using amp info for %s'%sensorId)
        detectors = self.camera.get(sensorId)

        # iterate over amplifiers
        for amp in detectors:
            ampName = amp.getName()[1:]
            hdu = segments[ampName]
            hdu.header['EXTNAME'] = 'Segment%s' % ampName
            hdu.header['DATASEC'] = noao_section_keyword(amp.getRawDataBBox())
            hdu.header['DETSEC'] = noao_section_keyword(
                amp.getBBox(),
                flipx=amp.getRawFlipX(),
                flipy=amp.getRawFlipY())

            # Flip up-down the intra-images for correct orientation
            if ampName in ['00', '01', '02', '03', '04', '05', '06', '07']:
                if verbose:
                    print('Flipping left-right..')
                data = np.copy(hdu.data)
                hdu.data[:] = np.fliplr(data)

            # Remove the incorrect BIASSEC keyword that phosim writes.
            try:
                hdu.header.remove('BIASSEC')
            except KeyError:
                pass
            sensor.append(hdu)

        # Set keywords in primary HDU, extracting most of the relevant
        # ones from the first phosim amplifier file.

        sensor[0].header['EXPTIME'] = sensor[1].header['EXPTIME']
        sensor[0].header['DARKTIME'] = sensor[1].header['DARKTIME']
        sensor[0].header['RUNNUM'] = sensor[1].header['OBSID']
        sensor[0].header['MJD-OBS'] = sensor[1].header['MJD-OBS']
        DATEOBS = astropy.time.Time(sensor[1].header['MJD-OBS'], format='mjd').isot
        sensor[0].header['DATE-OBS'] = DATEOBS
        YEAR, MONTH, DAYTIME = DATEOBS.split('-')
        DAY = DAYTIME[:2]
        DAYOBS = '%s%s%s'%(YEAR, MONTH, DAY)
        sensor[0].header['DAYOBS'] = DAYOBS

        DATE = sensor[1].header['DATE']  # file creation date
        sensor[0].header['DATE'] = DATE
        sensor[0].header['MJD'] = astropy.time.Time(DATE, format='isot').mjd

        sensor[0].header['FILTER'] = sensor[1].header['FILTER']

        serial = detectors.getSerial()  # eg. ITL-4400B-029
        CCD_MANU, CCD_TYPE, CCD_NUM = serial.split('-')
        sensor[0].header['LSST_NUM'] = serial
        sensor[0].header['CCD_MANU'] = CCD_MANU  # eg. ITL
        sensor[0].header['CCD_TYPE'] = CCD_TYPE  # eg. 4400B

        sensor[0].header['DETSIZE'] = noao_section_keyword(detectors.getBBox())
        sensor[0].header['INSTRUME'] = 'lsstCam'
        sensor[0].header['TELESCOP'] = 'LSST'

        sensor[0].header['TELCODE'] = TELCODE
        sensor[0].header['CONTRLLR'] = CONTRLLR
        SEQNUM = int(sensor[0].header['OBSID'][-6:])
        sensor[0].header['SEQNUM'] = SEQNUM

        OBSID = "%s_%s_%s_%s"%(TELCODE, CONTRLLR, DAYOBS, str(SEQNUM).zfill(6))
        sensor[0].header['OBSID'] = OBSID

        sensor[0].header['TESTTYPE'] = 'PHOSIM'
        sensor[0].header['IMGTYPE'] = 'SKYEXP'

        sensor[0].header['RAFTBAY'] = raft
        sensor[0].header['CCDSLOT'] = ccdslot

        sensor[0].header['RASTART'] = sensor[1].header['RA_DEG']
        sensor[0].header['DECSTART'] = sensor[1].header['DEC_DEG']
        sensor[0].header['ROTSTART'] = sensor[1].header['ROTANG']
        sensor[0].header['RA'] = sensor[1].header['RA_DEG']
        sensor[0].header['DEC'] = sensor[1].header['DEC_DEG']
        sensor[0].header['ROTPA'] = sensor[1].header['ROTANG']
        sensor[0].header['ROTPOS'] = sensor[1].header['ROTANG']

        # Filename created from TELCODE, CONTRLLR, DAYOBS, SEQNUM , raft, ccdslot.
        # eg. MC_C_20200825_000032_R00_SW0.fits
        # Below, OBSID already contains TELCODE, CONTRLLR, DAYOBS, str(SEQNUM).zfill(6)
        filename = '%s_%s_%s.fits'%(OBSID, raft, ccdslot)
        filename = os.path.join(out_dir, filename)

        # save the FITS file
        sensor.writeto(filename, overwrite=True)
        print('Saved as %s'%filename)

    def process_visit_eimage(self, visit_dir, out_dir=None, prefix='lsst',
                             verbose=False):
        """Process the visit eimage data.

        Parameters
        ----------
        visit_dir: str
            Directory containing the phosim eimage for a given visit.
        out_dir : str, optional
            Output directory for repackaged files. If None, then a directory
            with name 'eimage' will be created in the cwd. (the default is
            None.)
        prefix : str, optional
            Prefix to use to find PhoSim eimage files in visit_dir. (the
            default is 'lsst'.)
        verbose : bool, optional
            Set to True to print out time for processing each sensor. (the
            default is False.)
        """

        phosim_eimg_files \
            = sorted(glob.glob(os.path.join(visit_dir, f'{prefix}_e_*')))

        if (out_dir is None):
            out_dir = 'eimage'

        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        for eimg_file in phosim_eimg_files:
            if verbose:
                sys.stdout.write(eimg_file + '  ')
            t0 = time.time()

            print("repackaging", eimg_file)
            self.repackage_eimage(eimg_file, out_dir=out_dir)
            if (verbose):
                print(time.time() - t0)
                sys.stdout.flush()

    def repackage_eimage(self, phosim_eimg_file, out_dir='.', verbose=False):
        """Repackage the phosim eimage file to the format required by obs_lsst.

        Parameters
        ----------
        phosim_eimg_file : str
            PhoSim eimage file.
        out_dir : str, optional
            Output directory (the default is '.'.)
        verbose: bool, optional
            Set to True to print out processing information.
        """
        TELCODE = self.telcode
        CONTRLLR = self.contrllr

        # For corner rafts need to switch to
        # obs_lsst names for half-sensors
        ccd_mapper = {'C0': 'SW1',
                      'C1': 'SW0'}

        sensor = fits.open(phosim_eimg_file)[0]

        filename = phosim_eimg_file
        raft = os.path.basename(filename).split('_')[4]
        ccdslot = os.path.basename(filename).split('_')[5]
        detName = '%s_%s'%(raft, ccdslot)

        # Change the half-sensor name for corner rafts
        if detName in ['R00_S22', 'R04_S20', 'R40_S02', 'R44_S00']:
            ccdslot = os.path.basename(filename).split('_')[6]
            ccdslot = ccd_mapper[ccdslot]
            print('\nRepackaging corner ', raft, ccdslot)
        else:
            print('\nRepackaging main ', raft, ccdslot)

        # add required headers
        # eg. R00_SW1 (corner )  or R22_S00  (main)
        sensorId = '%s_%s'%(raft, ccdslot)

        detectors = self.camera.get(sensorId)

        # Set BOT-like keywords in the HDU
        sensor.header['EXPTIME'] = sensor.header['EXPTIME']
        sensor.header['DARKTIME'] = sensor.header['DARKTIME']
        sensor.header['RUNNUM'] = sensor.header['OBSID']
        sensor.header['MJD-OBS'] = sensor.header['MJD-OBS']

        DATEOBS = astropy.time.Time(sensor.header['MJD-OBS'], format='mjd').isot
        sensor.header['DATE-OBS'] = DATEOBS

        YEAR, MONTH, DAYTIME = DATEOBS.split('-')
        DAY = DAYTIME[:2]
        DAYOBS = '%s%s%s'%(YEAR, MONTH, DAY)
        sensor.header['DAYOBS'] = DAYOBS

        # file creation date as mjd
        sensor.header['MJD'] = astropy.time.Time(sensor.header['DATE'],
                                                 format='isot').mjd

        serial = detectors.getSerial()  # eg. ITL-4400B-029
        CCD_MANU, CCD_TYPE, CCD_NUM = serial.split('-')
        sensor.header['LSST_NUM'] = serial
        sensor.header['CCD_MANU'] = CCD_MANU  # eg. ITL
        sensor.header['CCD_TYPE'] = CCD_TYPE  # eg. 4400B

        sensor.header['DETSIZE'] = noao_section_keyword(detectors.getBBox())
        sensor.header['INSTRUME'] = 'lsstCam'
        sensor.header['TELESCOP'] = 'LSST'

        sensor.header['TELCODE'] = TELCODE
        sensor.header['CONTRLLR'] = CONTRLLR
        SEQNUM = int(sensor.header['OBSID'][-6:])
        sensor.header['SEQNUM'] = SEQNUM

        OBSID = "%s_%s_%s_%s"%(TELCODE, CONTRLLR, DAYOBS, str(SEQNUM).zfill(6))
        sensor.header['OBSID'] = OBSID

        sensor.header['TESTTYPE'] = 'PHOSIM'
        sensor.header['IMGTYPE'] = 'SKYEXP'

        sensor.header['RAFTBAY'] = raft
        sensor.header['CCDSLOT'] = ccdslot

        sensor.header['RASTART'] = sensor.header['RA_DEG']
        sensor.header['DECSTART'] = sensor.header['DEC_DEG']
        sensor.header['ROTSTART'] = sensor.header['ROTANG']

        sensor.header['RA'] = sensor.header['RA_DEG']
        sensor.header['DEC'] = sensor.header['DEC_DEG']
        sensor.header['ROTPA'] = sensor.header['ROTANG']
        sensor.header['ROTPOS'] = sensor.header['ROTANG']

        # Transpose the image to fulfill the geometry of postISRCCD by obs_lsst
        sensor.data = sensor.data.T

        # get a filename from TELCODE, CONTRLLR,  DAYOBS, SEQNUM , raft, ccdslot
        # eg. MC_H_20200825_000032_R00_SW0.fits
        filename = '%s_%s_%s.fits'%(OBSID, raft, ccdslot)
        filename = os.path.join(out_dir, filename)

        # save the FITS file
        sensor.writeto(filename, overwrite=True)
        print('Saved as %s'%filename)
