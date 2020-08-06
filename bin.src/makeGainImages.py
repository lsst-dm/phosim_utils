#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011, 2012, 2013 LSST Corporation.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
from datetime import datetime
import argparse

from lsst.obs.lsst.phosim import PhosimMapper
import lsst.afw.image as afwImage
import numpy as np

def main(just_wfs=False, detector_list=None):

    if (just_wfs is True) and (detector_list is not None):
        raise RuntimeError("--just_wfs and --detector_list are exclusive.")

    camera = PhosimMapper().camera
    if just_wfs:
        ccd_list = [camera[name] for name in ["R00_S22", "R04_S20", "R44_S00", "R40_S02"]]
    elif (detector_list is not None):
        ccd_list = [camera[name] for name in detector_list]
    else:
        ccd_list = camera

    for filt_name in 'ugrizy':
        for ccd in ccd_list:
            name = ccd.getName()
            # I'm not sure how to deal with the split chips yet.
            if 'A' in name or 'B' in name:
                continue
            print(name)
            CHIPID = "".join([c for c in name if c != "," and c != ":"])
            CHIPID = "_".join(CHIPID.split())
            image = afwImage.ImageF(ccd.getBBox())
            for amp in ccd:
                subim = afwImage.ImageF(image, amp.getBBox())
                subim[:] = 1/amp.getGain()
                print(amp.getName(), amp.getGain())
            
            # need to flip the image to match the result of phosim repackager... 
            oldImageArray = image.array.copy()
            image.array[:] = np.flipud(oldImageArray)

            expInfo = afwImage.ExposureInfo()
            inFilter = afwImage.Filter(filt_name)
            expInfo.setFilter(inFilter)
            exp = afwImage.ExposureF(afwImage.MaskedImageF(image), expInfo)
            md = exp.getMetadata()
            md.set('CHIPID', CHIPID)
            # Set place holder date
            md.set('MJD-OBS', 53005.0)
            md.set('OBSTYPE', 'flat')
            # arbitrary for flats
            md.set('EXPTIME', 100)
            # need to be able to specify any filter
            md.set('CALDATE', 53005.0)
            # Add the CALIB_ID header card
            md.set('CALIB_ID', 'raftName=%s detectorName=%s detector=%i filter=%s calibDate=%s' %
                   (CHIPID.split('_')[0], CHIPID.split('_')[1], ccd.getId(), filt_name, datetime.now()))
            exp.setMetadata(md)
            exp.writeFits("%(name)s_%(filter)s.fits"%({'name': CHIPID, 'filter': filt_name}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make fake flats using image gains')
    parser.add_argument('--just_wfs', action='store_true',
                        help='Generate fake flats for just wavefront sensing chips.')
    parser.add_argument('--detector_list', nargs='+',
                        help='Generate fake flats for the detector list. (e.g. R22_S11 R22_S10)')

    args = parser.parse_args()
    main(just_wfs=args.just_wfs, detector_list=args.detector_list)
